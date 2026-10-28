// Процесс SALMON_QC: сводная таблица проверок квантификации по всем образцам
// (mapping rate, фрагменты, тип библиотеки, пол, «домашние» гены, MT-риды + флаги).
// Файл *_mqc.tsv MultiQC подхватит как отдельный раздел отчёта.

process SALMON_QC {
    // в образе MultiQC уже есть Python — отдельный контейнер не нужен (linux/amd64)
    container 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'
    cpus 1

    input:
    // папки salmon всех образцов; stageAs кладёт их в подпапку salmon/ рабочей директории:
    //   salmon/N3/, salmon/T3/, ... — как ждёт скрипт
    path quant_dirs, stageAs: 'salmon/*'
    path tx2gene                     // tx2gene_gencodeNN.tsv

    output:
    path "salmon_qc_mqc.tsv", emit: tsv

    script:
    // salmon_qc.py лежит в bin/ проекта — Nextflow сам добавляет эту папку в PATH
    """
    salmon_qc.py salmon ${tx2gene} salmon_qc_mqc.tsv
    """
}