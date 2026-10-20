// Подключаем модули. FastQC — дважды под разными именами: до и после обрезки
include { FASTQC as FASTQC_RAW  } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { FASTP                 } from './modules/fastp'
include { MULTIQC               } from './modules/multiqc'

workflow {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists: true)

    FASTQC_RAW(reads_ch)                 // качество сырых ридов
    FASTP(reads_ch)                      // обрезка
    FASTQC_TRIM(FASTP.out.reads)         // качество после обрезки

    // всё для MultiQC — в один канал: mix сливает каналы, collect ждёт все задачи
    multiqc_files = FASTQC_RAW.out.zip
        .mix(FASTP.out.json, FASTQC_TRIM.out.zip)
        .collect()

    MULTIQC(multiqc_files)               // один общий отчёт
}