// Процесс FASTQC: отчёт о качестве для одной пары ридов
process FASTQC {
    tag "${sample}"
    container 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 2

    input:
    tuple val(sample), path(reads)

    output:
    path "*_fastqc.html", emit: html   // emit — имя выхода, чтобы обращаться к нему из workflow
    path "*_fastqc.zip",  emit: zip    // zip-архивы понадобятся MultiQC

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}

// Процесс MULTIQC: один сводный отчёт по всем образцам
process MULTIQC {
    container 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files              // список с путями к файлам

    output:
    path "multiqc_report.html" // сам отчёт
    path "multiqc_data"        // таблицы с цифрами из отчёта

    script:
    """
    multiqc .
    """
}

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC(reads_ch)

    // FASTQC выдаёт zip-архивы по одному на задачу; collect() ждёт,
    // пока закончатся ВСЕ задачи, и собирает архивы в один список
    MULTIQC(FASTQC.out.zip.collect())
}