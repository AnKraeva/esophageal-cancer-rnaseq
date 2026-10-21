// Подключаем модули. FastQC — дважды под разными именами: до и после обрезки
include { DOWNLOAD_ENA          } from './modules/download_ena'
include { FASTQC as FASTQC_RAW  } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { FASTP                 } from './modules/fastp'
include { MULTIQC               } from './modules/multiqc'

workflow {

    // samplesheet -> по строке на образец: [имя, ссылка R1, ссылка R2]
    samples_ch = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.fastq_1, row.fastq_2) }

    DOWNLOAD_ENA(samples_ch)             // риды из ENA (N пар или целиком)

    // [имя, R1, R2] -> [имя, [R1, R2]]: в таком виде риды ждут FastQC и fastp
    reads_ch = DOWNLOAD_ENA.out.reads
        .map { sample, r1, r2 -> tuple(sample, [r1, r2]) }

    FASTQC_RAW(reads_ch)                 // качество сырых ридов
    FASTP(reads_ch)                      // обрезка
    FASTQC_TRIM(FASTP.out.reads)         // качество после обрезки

    // всё для MultiQC — в один канал: mix сливает каналы, collect ждёт все задачи
    multiqc_files = FASTQC_RAW.out.zip
        .mix(FASTP.out.json, FASTQC_TRIM.out.zip)
        .collect()

    MULTIQC(multiqc_files)               // один общий отчёт
}