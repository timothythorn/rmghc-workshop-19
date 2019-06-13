#!/usr/bin/env nextflow
params.name             = "RNA-seq"
params.email            = "timothy.thorn@colorado.edu"
params.reads            = "/home/timothythorn/workshop/fastq/*{*_R1,*_R2}.fastq.gz"



log.info "RNA-seq Pipeline"
log.info "====================================="
log.info "name         : ${params.name}"
log.info "email        : ${params.email}"
log.info "reads        : ${params.reads}"
log.info "\n"



annotation = Channel.fromPath("/data/annotation/*")
genome = Channel.fromPath("/data/genome/*")
sample_info = Channel.fromPath("/data/sample_info/*")
index = Channel.fromPath("/data/index")
reads = Channel.fromFilePairs(params.reads, size: -1)
  .ifEmpty { error "Can't find any reads matching: ${params.reads}" }
  .into {
    reads_for_fastqc;
    reads_for_mapping
  }



annotation.into {
  annotation_for_count;
  annotation_for_transcriptome;
  annotation_for_de;
}



process make_transcriptome {

  maxForks 1
  publishDir 'results/genome'

  input:
  file annotation from annotation_for_transcriptome
  file genome from genome
  file sample_info from sample_info

  output:
  file "transcriptome.fa" into transcriptome
  file "gencode.vM21.annotation.gtf"
  file "sample_sheet.csv"

  script:
  """
  gffread -w transcriptome.fa -g ${genome} ${annotation}
  """
}



process fastqc {

  maxForks 1
  publishDir 'results/fastqc'

  input:
  set sample_id, file(fastqz) from reads_for_fastqc

  output:
  file "*.zip" into fastqc

  script:
  """
  fastqc --threads 4 -f fastq -q ${fastqz}
  """
}



process map {

  maxForks 1
  publishDir 'results/bam'

  input:
  set sample_id, file(reads), file(index) from reads_for_mapping.combine(index)

  output:
  set sample_id, file("*Aligned.out.bam") into mapped_genome
  set sample_id, file("*toTranscriptome.out.bam") into mapped_transcriptome
  file '*' into star

  script:
  """
  STAR  --runThreadN 5 \
  --genomeDir ${index} \
  --readFilesIn ${reads.findAll{ it =~ /\_R1\./ }.join(',')} \
                ${reads.findAll{ it =~ /\_R2\./ }.join(',')} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outSAMattributes NH HI NM MD AS \
  --outReadsUnmapped Fastx \
  --quantMode TranscriptomeSAM \
  --outFileNamePrefix ${sample_id}_
  """
}



mapped_genome.into {
  mapped_for_count;
  mapped_for_igv;
  mapped_for_markduplicates
}



process count {

  publishDir 'results/feature_counts'

  input:
  file annotation from annotation_for_count
  set sample_id, file(bam) from mapped_for_count

  output:
  file '*.fCounts'
  file '*.fCounts*' into counts

  script:
  """
  featureCounts  -C \
    -p \
    -T 1 \
    -g gene_id \
    -a ${annotation} \
    -o ${sample_id}.fCounts \
    ${bam}
  """
}



process salmon {

  publishDir 'results/salmon'

  input:
  file transcript_fasta from transcriptome
  set sample_id, file(bam) from mapped_transcriptome

  output:
  file("*") into salmon
  file("*") into salmon_for_de

  script:
  """
  salmon quant -l A \
    -p 1 \
    -t ${transcript_fasta} \
    -o ${sample_id} \
    -a ${bam} \
    --numBootstraps 30
  """
}



process sort_bam {

  publishDir 'results/igv'

  input:
  set sample_id, file(bam_file) from mapped_for_igv

  output:
  set sample_id, file("*.bam"), file('*.bai') into sorted_bam

  script:
  """
  samtools sort --threads 1 \
    -m 4G \
    -o ${sample_id}.bam \
    ${bam_file}
  samtools index ${sample_id}.bam
  """
}



process compile_fastqc {

  input:
  file fastqc from fastqc.collect()

  output:
  file "fastqc" into fastqc_compiled

  script:
  """
  mkdir fastqc
  mv ${fastqc} fastqc/.
  """
}



process compile_star {

  input:
  file star from star.collect()

  output:
  file "star" into star_compiled

  script:
  """
  mkdir star
  mv ${star} star/.
  """
}



process compile_salmon {

  input:
  file salmon from salmon.collect()

  output:
  file "salmon" into salmon_compiled

  script:
  """
  mkdir salmon
  mv ${salmon} salmon/.
  """
}



process compile_counts {

  input:
  file counts from counts.collect()

  output:
  file "counts" into counts_compiled

  script:
  """
  mkdir counts
  mv ${counts} counts/.
  """
}



process multiqc {

  publishDir "results/multiqc"

  input:
  file fastqc from fastqc_compiled
  file star from star_compiled
  file salmon from salmon_compiled
  file counts from counts_compiled

  output:
  set file('*_multiqc_report.html'), file('*_data/*')

  script:
  """
  multiqc ${fastqc} \
    ${star} \
    ${salmon} \
    ${counts} \
    --title '${params.name}' \
    --cl_config "extra_fn_clean_exts: [ '_1', '_2' ]"
  """
}



// process differential_expression {
//
//   publishDir 'reports', mode: 'copy'
//
//   input:
//   file sample_file from salmon_for_de.collect()
//   file annotation from annotation_for_de
//
//   script:
//   """
//   Rscript -e 'rmarkdown::render("${baseDir}/bin/differential_expression.Rmd", \
//     params = list(baseDir = "${baseDir}",\
//                   annotation_file = "${annotation}"))'
//   """
// }
