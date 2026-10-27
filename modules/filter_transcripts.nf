// Процесс FILTER_TRANSCRIPTS: оставляет в transcripts.fa только транскрипты основных хромосом
// (те, что есть в GTF). Убирает копии генов на альтернативных гаплотипах и патчах,
// иначе EM делит риды между копиями, а tximport потом выбрасывает «альтернативные»

process FILTER_TRANSCRIPTS {
    tag "GENCODE v${version}"
    // нужны только gzip и awk
    container 'community.wave.seqera.io/library/curl_gzip:57ce4c5aaf0198a6'
    cpus 1

    input:
    path transcripts                 // gencode.vNN.transcripts.fa.gz (все участки генома)
    path tx2gene                     // tx2gene_gencodeNN.tsv (только основные хромосомы)
    val version

    output:
    path "gencode.v${version}.transcripts.chr.fa.gz", emit: fasta

    script:
    // awk читает два файла по очереди:
    //   1-й (tx2gene): запоминаем ID транскриптов (пропуская строку-заголовок)
    //   2-й (FASTA, "-" = поток из gzip): на строке ">ENST...|ENSG...|..." берём ID до первого "|"
    //       и решаем, печатать ли эту запись (заголовок + строки последовательности под ним)
    """
    gzip -dc ${transcripts} \\
        | awk 'NR == FNR { if (FNR > 1) keep[\$1]; next }
               /^>/      { split(substr(\$0, 2), h, "|"); p = (h[1] in keep) }
               p' ${tx2gene} - \\
        | gzip -c > gencode.v${version}.transcripts.chr.fa.gz
    """
}