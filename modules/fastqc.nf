// Процесс FASTQC: отчёт о качестве для одной пары ридов
process FASTQC {
    tag "${sample}"
    container 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4'
    cpus 2

    input:
    tuple val(sample), path(reads)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip    // zip-архивы понадобятся MultiQC

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}