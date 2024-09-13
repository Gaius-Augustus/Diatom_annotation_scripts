nextflow.enable.dsl=2

process agat {
  publishDir "$params.outdir", mode: 'copy'

  container '/isg/shared/databases/nfx_singularity_cache/depot.galaxyproject.org-singularity-agat-1.4.0--pl5321hdfd78af_0.img'

  memory '10 GB'
  cpus 4 
  executor 'slurm'
  clusterOptions '--qos=general'
  queue 'general'

  input:
  tuple val(id), path(gff)

  output:
  path('*.gff'), emit: annotation
   
  script:

  """
  agat_convert_sp_gxf2gxf.pl -g ${gff} -o ${id}.gff
  """
}

workflow {
  ch_gff = Channel.fromPath(params.gff).flatten().map { file -> tuple(file.baseName, file) }
  agat     ( ch_gff )
}
