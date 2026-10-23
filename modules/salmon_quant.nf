// Процесс SALMON_QUANT: считает, сколько фрагментов пришлось на каждый транскрипт
// Вход: очищенные риды образца (после fastp) + готовый индекс

process SALMON_QUANT {
    tag "${sample}"
    // salmon 2.8.0, linux/amd64 (тот же образ, что в SALMON_INDEX)
    container 'community.wave.seqera.io/library/salmon:2.8.0--37bedaf2df54f35d'
    cpus 4                           // 2 образца параллельно × 4 ядра = 8 (лимит executor)
    memory '6 GB'                    // индекс в памяти (~2 ГБ) + запас; 2 × 6 = 12 ГБ < 16

    input:
    tuple val(sample), path(reads)   // [имя, [R1, R2]] из FASTP.out.reads
    path index                       // папка индекса (канал-значение: один на все образцы)

    output:
    // вся папка образца: quant.sf + логи и json (их читает MultiQC)
    tuple val(sample), path("${sample}"), emit: results

    script:
    // -l          тип библиотеки: A — salmon определит сам (ждём ISF)
    // --gcBias    поправка на GC-состав фрагментов
    // --seqBias   поправка на неслучайные начала фрагментов
    """
    salmon quant \\
        -i ${index} \\
        -l ${params.lib_type} \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --gcBias --seqBias \\
        -p ${task.cpus} \\
        -o ${sample}
    """
}