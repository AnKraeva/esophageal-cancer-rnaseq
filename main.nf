// Подключаем модули. FastQC — дважды под разными именами: до и после обрезки
include { DOWNLOAD_ENA          } from './modules/download_ena'
include { DOWNLOAD_REFERENCE    } from './modules/download_reference'
include { SALMON_INDEX          } from './modules/salmon_index'
include { FASTQC as FASTQC_RAW  } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { FASTP                 } from './modules/fastp'
include { MULTIQC               } from './modules/multiqc'
include { SALMON_QUANT          } from './modules/salmon_quant'
include { TX2GENE               } from './modules/tx2gene'
include { FILTER_TRANSCRIPTS    } from './modules/filter_transcripts'
include { SALMON_QC             } from './modules/salmon_qc'


workflow {

    // --- референс: не зависит от образцов, готовится один раз ---
    // на вход простое значение (а не канал из файла) → выходы станут каналами-значениями:
    // их можно читать сколько угодно раз (индекс понадобится всем 16 образцам)
    DOWNLOAD_REFERENCE(params.gencode_version)
    TX2GENE(DOWNLOAD_REFERENCE.out.gtf, params.gencode_version)   // транскрипт → ген (только основные хромосомы)
    // транскрипты только основных хромосом — ровно тот же набор, что в tx2gene
    FILTER_TRANSCRIPTS(DOWNLOAD_REFERENCE.out.transcripts, TX2GENE.out.tsv, params.gencode_version)
    SALMON_INDEX(FILTER_TRANSCRIPTS.out.fasta, params.gencode_version)

    // --- образцы ---
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




    // квантификация: каждый образец + один и тот же индекс
    SALMON_QUANT(FASTP.out.reads, SALMON_INDEX.out.index)

        // сводка проверок по всем образцам: собираем папки salmon в один список (collect)
    SALMON_QC(
        SALMON_QUANT.out.results.map { sample, dir -> dir }.collect(),
        TX2GENE.out.tsv
    )

    // всё для MultiQC — в один канал: mix сливает каналы, collect ждёт все задачи
    // из salmon берём только папку (второй элемент [имя, папка]) — MultiQC сам найдёт в ней логи
    multiqc_files = FASTQC_RAW.out.zip
        .mix(FASTP.out.json,
             FASTQC_TRIM.out.zip,
             SALMON_QUANT.out.results.map { sample, dir -> dir },
             SALMON_QC.out.tsv)
        .collect()

    MULTIQC(multiqc_files)               // один общий отчёт
}