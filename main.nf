// Подключаем модули под именами "для сырых ридов"
include { FASTQC  as FASTQC_RAW  } from './modules/fastqc'
include { MULTIQC as MULTIQC_RAW } from './modules/multiqc'

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC_RAW(reads_ch)

    // collect() ждёт все задачи FastQC и собирает архивы в один список
    MULTIQC_RAW(FASTQC_RAW.out.zip.collect())
}