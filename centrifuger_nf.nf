#!/usr/bin/env nextflow

/*Author:
* - Kiran Adhikari */

//***parameters***//
/* Default Pipeline parameters. They can be overridden using command line*/

//reads is input fastq file
params.reads = "/home/kiran/kiran/nextflow/nf-training/centrifuger/reads/ont_shotgun_combined.fastq"

// index is chertrifuger ref indexs database 
params.indexs = "/home/kiran/kiran/nextflow/nf-training/centrifuger/cfr_indexes"

// To save output in results dir 
params.outdir = "results"

//println "reads: $params.reads"

//*** Session Info***//

log.info """\
          METAGENOMICS NEXTFLOW PIPELINE
          ==============================
          indexs  : ${params.indexs}
          reads   : ${params.reads}
          outdir  : ${params.outdir}
          """
          .stripIndent()
          
//*** Channels***//

// reads are input fastq file, 
reads_ch = Channel.fromPath(params.reads, checkIfExists:true)
reads_ch.view()

//index are centrifuger reference indexes
indexs_ch = Channel.fromPath(params.indexs, checkIfExists:true)
indexs_ch.view()


//*** workflow ***//

workflow {

  main:
  centri_ch = METANALYSIS(indexs_ch, reads_ch)
  kreport_ch = KREPORT(centri_ch)
  
}

//***Process***//
//***** RUN centrifuger tool***//

process METANALYSIS {
  
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
  path indexs
  path reads
  
  output:
  path "reports/output.tsv"
  
  script:
  """
  mkdir -p reports
  centrifuger -x ${indexs}/cfr_ref_idx -u ${reads} -t 8 > reports/output.tsv
  """
}

//***Generate Kreport***///

process KREPORT {

  publishDir "${params.outdir}", mode: 'copy'
  
  input:
  path "reports/output.tsv"
  
  output:
  path "reports/kreport.centr.tsv"
  path "reports/final.kreport.tsv"
  
  script:
  """
  centrifuger-kreport -x ${params.indexs}/cfr_ref_idx reports/output.tsv > reports/kreport.centr.tsv

  # Add header to the kreport
  echo -e "%_clade_rooted_reads\\tNo._reads_assign_to_root\\tno.ofreads_assign_to_taxon\\tRank\\tTax_ID\\tScientific_Name" > reports/final.kreport.tsv

  # Filter reads from kreport to final report
  awk -F'\\t' '\$1 >=0.5 && \$4 == "S"' reports/kreport.centr.tsv  >> reports/final.kreport.tsv
  """
}




