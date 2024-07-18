#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params: parameters for workflow 
params.sleep = 2 
params.input = "/home/kiran/kiran/nextflow/nf-training/centrifuger/read_pod5"
params.model = "sup@v4.2.0"

// channels : to send data to workflow 
input_ch = Channel.fromPath(params.input)

//workflow 

workflow {
  
  BASECALLING(input_ch)
  
  BASECALLING.out.view()

 }

//process

// input just from Path
// output name of the output file 

process BASECALLING {
  
  input:
  path pod5
  
  output:
  path "basecall.fastq"
  
  script:
  """
  ### running dorado tool for basecalling 
  dorado basecaller ${params.model} ${pod5} --emit-fastq -r > basecall.fastq
  
  """
  
 }
 
 
 
 
 