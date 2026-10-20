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
    // префикс задаётся в конфиге (raw / trim); если не задан — просто имя образца
    def prefix = task.ext.prefix ?: sample
    """
    # ln -s создаёт ссылку на файл под новым именем (копировать не нужно)
    ln -s ${reads[0]} ${prefix}_R1.fastq.gz
    ln -s ${reads[1]} ${prefix}_R2.fastq.gz
    fastqc --threads ${task.cpus} ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz
    """
}