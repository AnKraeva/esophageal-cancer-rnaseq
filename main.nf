// Процесс FASTQC: отчёт о качестве для одной пары ридов
process FASTQC {
    tag "${sample}"   // подпись в логе: видно, какой образец сейчас обрабатывается
    container 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4'
    publishDir "${params.outdir}/fastqc", mode: 'copy'   // скопировать готовые отчёты в results/fastqc
    cpus 2

    input:
    tuple val(sample), path(reads)   // имя образца + список из двух файлов

    output:
    path "*_fastqc.html"   // отчёты для просмотра
    path "*_fastqc.zip"    // те же данные в машинном виде (понадобятся для MultiQC)

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}

// Главный workflow: откуда берём данные и куда их отправляем
workflow {
    // ищем файлы по шаблону и группируем по общей части имени:
    // SRR13344601 -> [SRR13344601_1.sub1M.fastq.gz, SRR13344601_2.sub1M.fastq.gz]
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)

    reads_ch.view()   // напечатать содержимое канала: удобно для отладки

    FASTQC(reads_ch)
}