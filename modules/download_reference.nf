// Процесс DOWNLOAD_REFERENCE: скачивает референс GENCODE для человека
//   transcripts.fa.gz — последовательности транскриптов (из них salmon строит индекс)
//   annotation.gtf.gz — аннотация (из неё сделаем tx2gene: транскрипт → ген)

process DOWNLOAD_REFERENCE {
    tag "GENCODE v${version}"
    // тот же контейнер, что в DOWNLOAD_ENA: curl + gzip, linux/amd64
    container 'community.wave.seqera.io/library/curl_gzip:57ce4c5aaf0198a6'
    cpus 1

    input:
    val version                      // номер релиза GENCODE, например 50

    output:
    // два отдельных выхода: транскрипты пойдут в индекс, GTF — в tx2gene
    path "gencode.v${version}.transcripts.fa.gz", emit: transcripts
    path "gencode.v${version}.annotation.gtf.gz", emit: gtf

    script:
    def base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version}"
    """
    for f in transcripts.fa.gz annotation.gtf.gz; do
        out="gencode.v${version}.\$f"

        # -f: при ошибке сервера (например, 404) не сохранять страницу ошибки как файл
        # --retry: файл качаем целиком, поэтому повтор безопасен (в отличие от DOWNLOAD_ENA)
        curl -fsSL --retry 5 --retry-delay 10 -o "\$out" "${base}/\$out"

        # проверка: архив целый, не оборван на середине
        gzip -t "\$out"
    done
    """
}