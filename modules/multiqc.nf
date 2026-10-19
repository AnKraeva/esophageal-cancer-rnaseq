// Процесс MULTIQC: один сводный отчёт по всем образцам
process MULTIQC {
    container 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'

    input:
    path qc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}