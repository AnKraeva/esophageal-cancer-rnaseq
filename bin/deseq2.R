# 1. Параметры запуска --------------------------------------------------------

# commandArgs() возвращает значения, которые передали после имени скрипта.
# Благодаря этим аргументам скрипт можно запустить из Jupyter или из Nextflow.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "Usage: Rscript bin/deseq2.R <samplesheet.csv> <tx2gene.tsv> ",
    "<salmon_dir> <output_dir> <gencode_version>",
    call. = FALSE
  )
}

samplesheet_path <- args[1]
tx2gene_path <- args[2]
salmon_dir <- args[3]
output_dir <- args[4]
gencode_version <- args[5]  # Записываем в run_info.txt.

# Параметры анализа из DESeq.ipynb.
min_count <- 10L
min_samples <- 8L
alpha <- 0.05
pca_n_genes <- 500L

# patient учитывает исходные различия между пациентами.
# condition оценивает эффект опухоли по сравнению с нормой.
design_formula <- ~ patient + condition
contrast <- c("condition", "tumor", "normal")
condition_coef <- "condition_tumor_vs_normal"  # Имя коэффициента для apeglm.

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})
if (!requireNamespace("apeglm", quietly = TRUE)) {
  stop("Install apeglm in the R environment before running this script.", call. = FALSE)
}

# Общие настройки для всех таблиц и графиков.
write_tsv <- function(data, filename) {
  write.table(
    data, file.path(output_dir, filename),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

save_png <- function(filename, draw, width = 8, height = 6) {
  # Открываем PNG, рисуем график и закрываем файл при выходе из функции.
  grDevices::png(
    file.path(output_dir, "plots", filename),
    width = width, height = height, units = "in", res = 160, type = "cairo"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
}


# 2. Образцы и входные файлы --------------------------------------------------

samples <- read.csv(samplesheet_path, colClasses = "character")

# По имени образца определяем пациента и ткань. Например, N3 и T3 –
# нормальная и опухолевая ткань пациента 3. Имена должны быть уникальными.
if (
  !"sample" %in% names(samples) ||
    anyNA(samples$sample) ||
    anyDuplicated(samples$sample) ||
    any(!grepl("^[NT][0-9]+$", samples$sample))
) {
  stop("The samplesheet must contain unique sample names like N3 and T3.", call. = FALSE)
}

samples$patient <- factor(sub("^[NT]", "", samples$sample))
samples$condition <- factor(
  ifelse(startsWith(samples$sample, "N"), "normal", "tumor"),
  levels = c("normal", "tumor")
)
rownames(samples) <- samples$sample

# Для парной модели у каждого пациента ждём ровно один normal и один tumor.
# Порог фильтрации также требует, чтобы всего было хотя бы 8 образцов.
if (nrow(samples) < min_samples || any(table(samples$patient, samples$condition) != 1L)) {
  stop(
    "DE requires at least ", min_samples,
    " samples with exactly one normal and one tumor sample per patient.",
    call. = FALSE
  )
}

# Имена вектора станут именами колонок после импорта tximport.
quant_files <- setNames(file.path(salmon_dir, samples$sample, "quant.sf"), samples$sample)
missing_quant <- !file.exists(quant_files)
if (any(missing_quant)) {
  stop(
    "Missing quant.sf for: ", paste(names(quant_files)[missing_quant], collapse = ", "),
    call. = FALSE
  )
}

# Для импорта нужны ID транскрипта и гена. Название гена добавим при экспорте.
annotation <- read.delim(tx2gene_path, colClasses = "character")
annotation_columns <- c("transcript_id", "gene_id", "gene_name")
if (!all(annotation_columns %in% names(annotation))) {
  stop("tx2gene must contain transcript_id, gene_id and gene_name.", call. = FALSE)
}

tx2gene <- annotation[, c("transcript_id", "gene_id")]
if (
  anyNA(tx2gene) ||
    any(!nzchar(trimws(tx2gene$transcript_id))) ||
    any(!nzchar(trimws(tx2gene$gene_id))) ||
    anyDuplicated(tx2gene$transcript_id)
) {
  stop("tx2gene contains missing IDs or duplicate transcript IDs.", call. = FALSE)
}


# 3. Импорт, фильтрация, нормализация и VST -----------------------------------

# tximport суммирует оценки числа фрагментов по транскриптам одного гена.
# Он также передаёт средние эффективные длины, которые нужны при нормализации.
txi <- tximport(
  quant_files, type = "salmon", tx2gene = tx2gene,
  countsFromAbundance = "no"
)

# Перед созданием объекта DESeq2 проверяем, что метаданные идут
# в том же порядке, что и колонки с каунтами.
stopifnot(identical(colnames(txi$counts), rownames(samples)))
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)
n_genes_before_filter <- nrow(dds)

# Для каждого гена считаем, в скольких образцах каунт >= 10.
# Оставляем гены, у которых таких образцов хотя бы 8.
keep <- rowSums(counts(dds) >= min_count) >= min_samples
if (!any(keep)) {
  stop("No genes passed the count filter.", call. = FALSE)
}

# DESeq2 рассчитает поправки на особенности библиотек и эффективные длины.
# VST использует эти поправки и даёт значения для отдельной задачи – PCA.
dds <- estimateSizeFactors(dds[keep, ])
vsd <- vst(dds, blind = FALSE)


# 4. Дифференциальная экспрессия и shrinkage ---------------------------------

# Модель DESeq2 получает отфильтрованный объект dds с каунтами и поправками.
# Тест Вальда даёт p-value для эффекта опухоли с учётом пациента.
dds_de <- DESeq(dds, test = "Wald")
res <- results(dds_de, contrast = contrast, alpha = alpha)

# apeglm уточняет величину эффекта для того же сравнения.
# Значимость по-прежнему берём из res, где лежат исходные pvalue и padj.
res_shrunk <- lfcShrink(dds_de, coef = condition_coef, res = res, type = "apeglm")

# Ставим обе оценки LFC в одну строку по ID гена. Символ гена добавляем
# из аннотации, а ID сохраняем с версией, как в исходных данных.
gene_ids <- rownames(res)
res_shrunk <- res_shrunk[gene_ids, ]
gene_names <- annotation$gene_name[match(gene_ids, annotation$gene_id)]

de_table <- data.frame(
  gene_id = gene_ids,
  gene_name = gene_names,
  baseMean = res$baseMean,
  log2FC_raw = res$log2FoldChange,
  lfcSE_raw = res$lfcSE,
  wald_stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  log2FC_apeglm = res_shrunk$log2FoldChange,
  lfcSD_apeglm = res_shrunk$lfcSE  # Для apeglm это апостериорное SD.
)

# padj определяет значимость, знак исходного LFC – направление,
# а модуль LFC после apeglm – порядок генов внутри каждого списка.
significant <- de_table[!is.na(de_table$padj) & de_table$padj < alpha, ]
significant <- significant[order(-abs(significant$log2FC_apeglm), significant$gene_id), ]
up <- subset(significant, log2FC_raw > 0)
down <- subset(significant, log2FC_raw < 0)


# 5. Экспорт таблиц ---------------------------------------------------------

# Из dds_de извлекаем каунты с нормировочными поправками.
# Строки идут как в de_results.tsv, колонки образцов – как в samplesheet.
normalized_counts <- counts(dds_de, normalized = TRUE)
counts_table <- data.frame(
  gene_id = gene_ids,
  gene_name = gene_names,
  normalized_counts[gene_ids, samples$sample, drop = FALSE],
  check.names = FALSE
)

dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
write_tsv(de_table, "de_results.tsv")
write_tsv(up, "significant_up.tsv")
write_tsv(down, "significant_down.tsv")
write_tsv(counts_table, "normalized_counts.tsv")


# 6. Графики ----------------------------------------------------------------

# PCA получает VST-значения. Каждая точка – один образец.
save_png("pca.png", function() {
  pca <- plotPCA(vsd, intgroup = "condition", ntop = pca_n_genes) +
    ggplot2::geom_text(
      ggplot2::aes(label = name), nudge_x = 3, vjust = -0.6,
      size = 3, show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 28, 8, 8))
  print(pca)
}, width = 10)

# Оценки дисперсии для генов и общая зависимость от средней экспрессии.
save_png("dispersion.png", function() {
  plotDispEsts(dds_de)
})

# MA: средняя экспрессия по X, LFC после apeglm по Y.
# Подсветка значимых генов использует исходный padj.
save_png("ma.png", function() {
  ma_res <- res_shrunk
  ma_res$padj <- res$padj
  plotMA(ma_res, alpha = alpha, main = "Tumor vs normal")
})

# Volcano: для координаты по Y берём -log10(padj). Если padj округлился
# до нуля, подставляем малое положительное число только для графика.
volcano <- subset(de_table, is.finite(log2FC_apeglm) & !is.na(padj))
volcano$y <- -log10(pmax(volcano$padj, .Machine$double.xmin))

# Подписываем до шести генов из каждого отсортированного списка.
# Для гена без названия используем его ID.
top_genes <- rbind(
  head(subset(up, is.finite(log2FC_apeglm)), 6L),
  head(subset(down, is.finite(log2FC_apeglm)), 6L)
)
labels <- top_genes$gene_name
missing_labels <- is.na(labels) | !nzchar(labels)
labels[missing_labels] <- top_genes$gene_id[missing_labels]

save_png("volcano.png", function() {
  if (nrow(volcano) == 0L) {
    plot.new()
    title(main = "Tumor vs normal", sub = "No genes with finite LFC and adjusted p-values")
    return(invisible(NULL))
  }
  plot(
    y ~ log2FC_apeglm, data = volcano,
    col = ifelse(volcano$padj < alpha, "blue", "grey70"),
    pch = 16, cex = 0.55, xlab = "log2 fold change (apeglm)",
    ylab = "-log10(padj)", main = "Tumor vs normal"
  )
  abline(h = -log10(alpha), lty = 2, col = "grey40")
  if (nrow(top_genes) > 0L) {
    text(
      top_genes$log2FC_apeglm, volcano$y[match(top_genes$gene_id, volcano$gene_id)],
      labels = labels, pos = ifelse(top_genes$log2FC_apeglm < 0, 4, 2),
      cex = 0.6, col = "grey45"
    )
  }
}, width = 10, height = 8)

# Красная диагональ показывает, где исходный и скорректированный LFC равны.
save_png("lfc_shrinkage.png", function() {
  plot(
    de_table$log2FC_raw, de_table$log2FC_apeglm,
    xlab = "Original log2 fold change", ylab = "log2 fold change (apeglm)",
    pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2)
  )
  abline(a = 0, b = 1, col = "red", lty = 2)
})


# 7. Параметры и окружение запуска ------------------------------------------

# Собираем в одном файле настройки анализа, числа для сверки с ноутбуком
# и версии R и пакетов, которые фактически использовались при запуске.
run_info <- c(
  "Analysis",
  paste("GENCODE release:", gencode_version),
  paste("Samples:", paste(samples$sample, collapse = ", ")),
  paste("Design:", deparse(design_formula)),
  paste("Contrast:", contrast[2], "vs", contrast[3]),
  "Import: Salmon, countsFromAbundance = no",
  sprintf("Count filter: >= %d counts in >= %d samples", min_count, min_samples),
  "VST: blind = FALSE",
  paste("PCA variable genes:", pca_n_genes),
  "Test: Wald; adjusted p-values: Benjamini-Hochberg with independent filtering",
  paste("Significance: padj <", alpha),
  "LFC shrinkage: apeglm",
  paste("Genes before count filter:", n_genes_before_filter),
  paste("Genes after count filter:", nrow(de_table)),
  paste("Significant genes:", nrow(significant)),
  paste("Higher in tumor:", nrow(up)),
  paste("Lower in tumor:", nrow(down)),
  "",
  "R session",
  capture.output(sessionInfo())
)

writeLines(run_info, file.path(output_dir, "run_info.txt"))
message("Saved DE outputs in ", output_dir)