nextflow.enable.dsl=2

process FETCH_SRA {
  publishDir "$params.outdir/sra", mode: 'copy'

  container '/isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-sra-tools-3.1.1--h4304569_2.img'

  memory '30 GB'
  cpus 8 

  input:
  val(sra)

  output:
  tuple val(sra), path("*.fastq"),  emit: sra
   
  script:

  """
  mkdir -p tmp
  export TMPDIR=tmp
  fasterq-dump ${sra} -e ${task.cpus}

  """
}
process KRAKEN {
publishDir "$params.outdir/kraken2", mode: 'copy'
  memory '70 GB'
  cpus 8 

container '/core/labs/Wegrzyn/diatom_genomes/smudgeplot/final_run/kraken-2.1.3.sif'
  input:
  tuple val(id), path(reads)
  
  output:
  tuple val(id), path("*unclassified*.fastq"), emit: reads
  tuple val(id), path("*.report"), emit: report
  tuple val(id), path("*output"), emit: output
  script:

"""
kraken2 --db /isg/shared/databases/kraken2/b+a+v+f --report ${id}.report --paired --threads 8 ${reads} --output ${id}.output --unclassified-out ${id}_unclassified#.fastq
"""
}
process FASTK {
  publishDir "$params.outdir/fastk", mode: 'copy'

  memory '16 GB'
  cpus 12 

  container '/isg/shared/databases/nfx_singularity_cache/smudgeplot-0.3.0.sif'

  input:
  tuple val(id), path(reads)
  val(kmer)

  output:
  tuple val(id), path("*${id}_FastK_Table*"), emit: ktab
  tuple val(id), path(".${id}_FastK_Table*"), emit: ktab_hidden

  script:
  """
  mkdir -p tmpdir
  export TMPDIR=tmpdir
  FastK -v -t4 -k${kmer} -M16 -T12 ${reads} -Ptmpdir -N"${id}"_FastK_Table

  """
}
process SMUDGEPLOT {
  publishDir "$params.outdir/smudgeplot", mode: 'copy'

  memory '30 GB'
  cpus 4 

  container '/isg/shared/databases/nfx_singularity_cache/smudgeplot-0.3.0.sif'

  input:
  tuple val(id), path(fastk)
  tuple val(id), path(hidden)

  output:
  path("${id}_run*"), emit: files

  script:
  """
  mkdir -p tmpdir
  export TMPDIR=tmpdir
  Histex -G ${id}_FastK_Table > ${id}.hist
  L=\$(smudgeplot.py cutoff ${id}.hist L)
  smudgeplot.py hetmers -L \$L -t 4 -o kmerpairs -tmp tmpdir --verbose ${id}_FastK_Table
  smudgeplot.py plot -o ${id}_run kmerpairs_text.smu

  """
}
workflow {
  sra_ch = Channel.fromPath(params.sra, checkIfExists: true).splitCsv().flatten()
  FETCH_SRA   ( sra_ch )
  KRAKEN     (FETCH_SRA.out.sra)
  FASTK       ( KRAKEN.out.reads, params.kmer )
  SMUDGEPLOT  ( FASTK.out.ktab, FASTK.out.ktab_hidden )
}
