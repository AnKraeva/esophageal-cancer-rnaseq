// Процесс TX2GENE: таблица «транскрипт → ген → имя гена» из GTF
// Нужна в R: tximport складывает транскрипты из quant.sf в гены

process TX2GENE {
    tag "GENCODE v${version}"
    // нужны только gzip и awk — берём тот же лёгкий образ, что для загрузки
    container 'community.wave.seqera.io/library/curl_gzip:57ce4c5aaf0198a6'
    cpus 1

    input:
    path gtf                         // gencode.vNN.annotation.gtf.gz
    val version

    output:
    path "tx2gene_gencode${version}.tsv", emit: tsv

    script:
    // берём строки с типом "transcript" (3-я колонка) и разбираем атрибуты (9-я колонка):
    //   gene_id "ENSG..."; transcript_id "ENST..."; gene_name "TP53"; ...
    // ID оставляем С ВЕРСИЯМИ (.1, .2...) — так же, как в quant.sf
    """
    printf "transcript_id\\tgene_id\\tgene_name\\n" > tx2gene_gencode${version}.tsv

    gzip -dc ${gtf} | awk -F'\\t' 'BEGIN { OFS = "\\t" }
        \$3 == "transcript" {
            tx = ""; gene = ""; name = ""
            n = split(\$9, attrs, ";")            # атрибуты разделены точкой с запятой
            for (i = 1; i <= n; i++) {
                sub(/^ +/, "", attrs[i])          # убрать пробел в начале
                split(attrs[i], kv, " ")          # kv[1] = ключ, kv[2] = "значение"
                v = kv[2]; gsub(/"/, "", v)       # убрать кавычки
                if      (kv[1] == "transcript_id") tx   = v
                else if (kv[1] == "gene_id")       gene = v
                else if (kv[1] == "gene_name")     name = v
            }
            print tx, gene, name
        }' >> tx2gene_gencode${version}.tsv
    """
}