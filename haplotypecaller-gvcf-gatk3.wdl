## Copyright Broad Institute, 2018
## 
## This WDL workflow runs HaplotypeCaller in VCF mode on a single sample according to 
## the GATK Best Practices (June 2016) scattered across intervals.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One VCF file and its index
##
## Cromwell version support 
## - Successfully tested on v33
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION 
workflow HaplotypeCallerGvcf_GATK3 {
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File scattered_calling_intervals_list
  
  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  Boolean? make_gvcf
  Boolean making_gvcf = select_first([make_gvcf, false])

  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename 
  String vcf_index = sample_basename 
  
  String output_suffix = if making_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_filename = vcf_basename + output_suffix

  #Docker
  String? gitc_docker_override
  String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])
  String? picard_docker_override
  String picard_docker = select_first([picard_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])
  String? picard_path_override
  String picard_path = select_first([picard_path_override, "/usr/gitc/"])

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = length(scattered_calling_intervals) - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval_list = interval_file,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_gvcf = making_gvcf,
        docker = gitc_docker,
    }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      docker = picard_docker,
      picard_path = picard_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }
}

# TASK DEFINITIONS

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  
  # Command Paramters
  File input_bam
  File input_bam_index
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File interval_list
  Int? contamination
  Boolean make_gvcf
  Int hc_scatter

  String? java_opt
  String java_option = select_first([java_opt,"-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m"])

  # Runtime parameters
  String docker
  Int? machine_mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  # We use PringReads to add interval_padding 500 to make sure that HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # This isn't needed with HaploypeCaller in  GATK4 because it can stream the required intervals directly from the cloud.

  command {

    /usr/gitc/gatk4/gatk-launch --javaOptions "-Xms2g" \
      PrintReads \
      -I ${input_bam} \
      --interval_padding 500 \
      -L ${interval_list} \
      -O local.sharded.bam \
    && \
    java ${java_option} \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${output_filename} \
      -I local.sharded.bam \
      -L ${interval_list} \
      --max_alternate_alleles 3 \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead \
      ${true="-ERC GVCF" false="" make_gvcf}

  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 7]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {

  # Command Paramters
  Array [File] input_vcfs
  Array [File] input_vcfs_indexes
  String output_filename
  Int compression_level
  String picard_path
  String java_opt

  # Runtime Paramters
  String docker
  Int? machine_mem_gb 
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int disk_size = ceil(size(write_lines(input_vcfs)) + 30) 

  command {
    java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${picard_path}picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_filename}
  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 3]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
}

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

