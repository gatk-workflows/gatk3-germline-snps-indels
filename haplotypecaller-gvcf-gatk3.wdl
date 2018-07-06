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
workflow HaplotypeCallerVcf_GATK3 {
  File input_bam
  File input_bam_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File scattered_calling_intervals_list
  
  Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

  Boolean? make_gvcf
  Boolean making_gvcf = select_first([make_gvcf,false])

  String sample_basename = basename(input_bam, ".bam")
  
  String vcf_basename = sample_basename 
  String vcf_index = sample_basename 
  
  String output_suffix = if making_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_filename = vcf_basename + output_suffix

  #Docker
  String? gatk_docker
  String gatk_image = select_first([gatk_docker, "broadinstitute/gatk3:3.8-1"])
  String? gatk_path
  String gatk_path2launch = select_first([gatk_path, "/usr/"])
  String? picard_docker
  String picard_image = select_first([picard_docker, "broadinstitute/genomes-in-the-cloud:2.3.0-1501082129"])
  String? picard_path
  String picard_path2launch = select_first([picard_path, "/usr/gitc/"])

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
        docker = gatk_image,
        gatk_path = gatk_path2launch
    }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indices = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      docker = picard_image,
      picard_path = picard_path2launch
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_merged_vcf = MergeVCFs.output_vcf
    File output_merged_vcf_index = MergeVCFs.output_vcf_index
  }
}

# TASK DEFINITIONS

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File interval_list
  Int? interval_padding
  Int? contamination
  Int? max_alt_alleles
  Boolean make_gvcf
  Int hc_scatter

  String gatk_path
  String java_opt

  # Runtime parameters
  String docker
  Int? mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts

  Int machine_mem_gb = select_first([mem_gb, 7])

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  command {
    java ${java_opt} -jar ${gatk_path}GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -o ${output_filename} \
      -L ${interval_list} \
      -ip ${default=100 interval_padding} \
      --max_alternate_alleles ${default=3 max_alt_alleles} \
      --read_filter OverclippedRead \
      -contamination ${default=0 contamination} \
      ${true="-ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000" false="" make_gvcf}
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
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
  Array [File] input_vcfs
  Array [File] input_vcfs_indices
  String output_filename

  Int compression_level
  Int? preemptible_attempts
  Int disk_size
  String mem_size

  String docker
  String picard_path
  String java_opt

  command {
    java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${picard_path}picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_filename}
  }

  runtime {
    docker: docker
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
}

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

