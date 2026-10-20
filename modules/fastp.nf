// Процесс FASTP: вырезает UMI и линкер, обрезает адаптеры, фильтрует риды
process FASTP {
    tag "${sample}"
    // контейнер на seqera.io/containers: поиск "fastp", платформа linux/amd64
    container 'community.wave.seqera.io/library/fastp:1.4.0--8808d8f6dbcbc281'
    cpus 4

    input:
    tuple val(sample), path(reads)

    output:
    // очищенная пара ридов вместе с именем образца – пойдёт в FastQC, потом в salmon
    tuple val(sample), path("${sample}_*.trim.fastq.gz"), emit: reads
    path "${sample}.json", emit: json   // для MultiQC
    path "${sample}.html", emit: html   // отдельный отчёт fastp


    script:
    // --umi ...                UMI: 12 нт в начале Read 1 → в имя рида, следом 13 нт линкера выбросить
    // --detect_adapter_for_pe  поиск адаптера по последовательности в дополнение к перекрытию пары

    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample}_1.trim.fastq.gz -O ${sample}_2.trim.fastq.gz \\
        --umi --umi_loc read1 --umi_len 12 --umi_skip 13 --umi_prefix UMI \\
        --detect_adapter_for_pe \\
        -w ${task.cpus} \\
        -j ${sample}.json -h ${sample}.html
    """
}