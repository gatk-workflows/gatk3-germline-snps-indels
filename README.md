# gatk3-germline-snps-indels

### Purpose : 
Workflows for germline short variant discovery with GATK3. 

### haplotypecaller-gvcf-gatk :
The haplotypecaller-vcf-gatk3 workflow runs HaplotypeCaller 
from GATK3 in VCF mode on a single sample according to the GATK Best Practices (June 2016), 
scattered across intervals.

#### Requirements/expectations
- One analysis-ready BAM file for a single sample (as identified in RG:SM)
- Set of variant calling intervals lists for the scatter, provided in a file
#### Outputs 
- One VCF file and its index

### Software version requirements :
- GATK 3 
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
- Successfully tested on v33
- Does not work on versions < v23 due to output syntax

### IMPORTANT NOTE : 
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
