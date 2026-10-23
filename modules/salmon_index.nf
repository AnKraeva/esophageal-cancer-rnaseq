// Процесс SALMON_INDEX: строит индекс salmon по транскриптам GENCODE (без декоев)
// Строится один раз, потом используется для всех образцов

process SALMON_INDEX {
    tag "GENCODE v${version}"
    container 'community.wave.seqera.io/library/salmon:2.8.0--37bedaf2df54f35d'
    cpus 8
    memory '8 GB'                    // для индекса только по транскриптам с запасом

    input:
    path transcripts                 // gencode.vNN.transcripts.fa.gz
    val version                      // номер GENCODE — для имени папки индекса

    output:
    // версия в имени: при переходе на новый GENCODE не подхватится старый индекс из storeDir
    path "salmon_index_gencode${version}", emit: index

    script:
    // -k 31      длина k-меров (значение по умолчанию, пишем явно для наглядности)
    // --gencode  обрезать заголовки GENCODE "ENST...|ENSG...|..." до ID транскрипта
    """
    salmon index \\
        -t ${transcripts} \\
        -i salmon_index_gencode${version} \\
        -k 31 \\
        --gencode \\
        -p ${task.cpus}
    """
}