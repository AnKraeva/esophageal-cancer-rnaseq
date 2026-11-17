process DESEQ2 {
    cpus 1
    container 'community.wave.seqera.io/library/bioconductor-apeglm_bioconductor-deseq2_bioconductor-tximport_jupyterlab_pruned:4135698481939e4f'
    memory '4 GB'

    input:
    // В рабочей папке получим salmon/N3/quant.sf, salmon/T3/quant.sf и т. д.
    path quant_dirs, stageAs: 'salmon/*'
    path tx2gene
    path samplesheet
    path r_script
    val gencode_version

    output:
    path 'de', emit: results

    script:
    """
    Rscript "${r_script}" "${samplesheet}" "${tx2gene}" salmon de "${gencode_version}"
    """
}