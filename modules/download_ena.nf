// Процесс DOWNLOAD_ENA: скачивает риды образца из ENA.
//   n_reads > 0 — только первые N пар ридов (читаем поток и останавливаемся на N)
//   n_reads = 0 — файлы целиком

// Число ридов из параметров. Из командной строки (--n_reads 100000) значение
// может прийти текстом "100000", поэтому явно превращаем его в число
def nReads() {
    return params.n_reads as long
}

// Метка для имени файла: 0 -> "full", 10000000 -> "sub10M", 100000 -> "sub100k"
def subLabel(long n) {
    if (n == 0)           return "full"
    if (n % 1000000 == 0) return "sub${n.intdiv(1000000)}M"
    if (n % 1000 == 0)    return "sub${n.intdiv(1000)}k"
    return "sub${n}"
}

process DOWNLOAD_ENA {
    tag "${sample}"
    // контейнер на seqera.io/containers: пакеты curl + gzip, платформа linux/amd64
    container 'community.wave.seqera.io/library/curl_gzip:57ce4c5aaf0198a6'
    cpus 2
    maxForks 4        // не больше 4 загрузок одновременно — бережём сеть и сервер ENA

    // Сервер ENA иногда обрывает соединение. Тогда задача перезапускается целиком
    // (до 3 раз) с паузой: 1 мин, 2 мин, 3 мин. Остальные загрузки при этом не прерываются
    errorStrategy { sleep(60000 * task.attempt); return 'retry' }
    maxRetries 3

    input:
    tuple val(sample), val(url_1), val(url_2)

    output:
    // метка в имени файла: другое N — другие файлы, старые не перезапишутся
    tuple val(sample),
          path("${sample}_1.${subLabel(nReads())}.fastq.gz"),
          path("${sample}_2.${subLabel(nReads())}.fastq.gz"), emit: reads

    script:
    def n_reads = nReads()
    def n_lines = n_reads * 4            // в FASTQ 4 строки на рид
    def label   = subLabel(n_reads)
    """
    for i in 1 2; do
        url=\$( [ \$i -eq 1 ] && echo "${url_1}" || echo "${url_2}" )
        out="${sample}_\${i}.${label}.fastq.gz"

        if [ ${n_reads} -eq 0 ]; then
            # --- весь файл: сохраняем .fastq.gz как есть ---
            curl -sSL --retry 5 --retry-delay 10 -o "\$out" "\$url"
            # проверка: архив целый, не оборван на середине
            gzip -t "\$out"
        else
            # --- первые N ридов ---
            # head закрывает поток, набрав нужное число строк; curl и gzip при этом
            # получают сигнал SIGPIPE и завершаются с ошибкой. Для этой цепочки это норма,
            # поэтому временно отключаем pipefail, а полноту файла проверяем ниже

            set +o pipefail
            # --retry здесь нельзя: curl начал бы файл заново и дописал бы его
            # в тот же поток → риды задвоились бы. Повтор делает Nextflow (errorStrategy)
            curl -sSL "\$url" \\
                | gzip -dc \\
                | head -n ${n_lines} \\
                | gzip -c > "\$out"
            set -o pipefail

            # проверка: в файле ровно N ридов (иначе загрузка оборвалась)
            got=\$(gzip -dc "\$out" | wc -l)
            if [ "\$got" -ne ${n_lines} ]; then
                echo "ОШИБКА: \$out — \$((got / 4)) ридов вместо ${n_reads}" >&2
                exit 1
            fi
        fi
    done
    """
}