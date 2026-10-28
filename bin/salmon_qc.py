#!/usr/bin/env python3
"""
Сводка проверок квантификации salmon по всем образцам — одна таблица.

Для каждого образца берёт:
  aux_info/meta_info.json   — сколько фрагментов обработано и картировано, длина фрагментов
  lib_format_counts.json    — тип библиотеки, совместимость, специфичность по цепи
  quant.sf + tx2gene        — транскрипты без гена, пол, «домашние» гены, доля MT-ридов
и помечает подозрительные значения в колонке flags.

Запуск (из папки проекта):
  python3 bin/salmon_qc.py results/salmon data/reference/tx2gene_gencode50.tsv results/qc/salmon_qc_mqc.tsv

Окончание _mqc.tsv и строки "# ..." в начале файла — чтобы MultiQC
подхватил таблицу как отдельный раздел отчёта.
Только стандартная библиотека Python, ничего устанавливать не нужно.
"""
import csv
import json
import statistics
import sys
from pathlib import Path

Y_GENES = ["RPS4Y1", "DDX3Y"]          # гены Y-хромосомы (мужской пол)
SEX_GENE_F = "XIST"                     # активен только у женщин
HOUSEKEEPING = ["ACTB", "GAPDH"]        # «домашние» гены: должны быть высоко


def load_tx2gene(path):
    """tx2gene.tsv → словарь {ID транскрипта: имя гена}."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {row["transcript_id"]: row["gene_name"] for row in reader}


def sum_by_gene(quant_path, tx2gene):
    """quant.sf → риды по генам, сумма всех ридов и число транскриптов без гена."""
    genes, total, no_gene = {}, 0.0, 0
    with open(quant_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            reads = float(row["NumReads"])
            gene = tx2gene.get(row["Name"])
            if gene is None:            # транскрипта нет в tx2gene → tximport его выбросит
                no_gene += 1
                continue
            genes[gene] = genes.get(gene, 0.0) + reads
            total += reads
    return genes, total, no_gene

def patient_order(name):
    """Ключ сортировки: сначала номер пациента (как число), внутри пары — N, потом T.
    N3, T3, N4, T4, ..., N22, T22. Имена не по схеме буква+число — в конец по алфавиту."""
    num = name[1:]
    return (0, int(num), name[0]) if num.isdigit() else (1, 0, name)

def call_sex(xist, y):
    """Пол по экспрессии: XIST — женский, гены Y — мужской. «?» если неясно."""
    if xist > 5 * (y + 1):
        return "F"
    if y > 5 * (xist + 1):
        return "M"
    return "?"


def collect(sample_dir, tx2gene):
    """Все метрики одного образца."""
    meta = json.loads((sample_dir / "aux_info" / "meta_info.json").read_text())
    lib = json.loads((sample_dir / "lib_format_counts.json").read_text())
    genes, total, no_gene = sum_by_gene(sample_dir / "quant.sf", tx2gene)

    xist = genes.get(SEX_GENE_F, 0)
    y = sum(genes.get(g, 0) for g in Y_GENES)
    mt = sum(v for g, v in genes.items() if g.startswith("MT-"))
    top_gene = max(genes, key=genes.get)

    return {
        "sample": sample_dir.name,
        "processed_M": meta["num_processed"] / 1e6,      # пары после fastp
        "mapped_M": meta["num_mapped"] / 1e6,            # картированные пары
        "mapping_pct": meta["percent_mapped"],
        "frag_mean": meta["frag_length_mean"],
        "lib_type": lib["expected_format"],
        "compat": lib["compatible_fragment_ratio"],
        "strand": lib["strand_mapping_bias"],
        "no_gene_tx": no_gene,
        "XIST": xist,
        "Y_genes": y,
        "sex": call_sex(xist, y),
        "ACTB": genes.get("ACTB", 0),
        "GAPDH": genes.get("GAPDH", 0),
        "MT_pct": 100 * mt / total if total else 0,
        "top_gene": top_gene,
        "top_gene_pct": 100 * genes[top_gene] / total if total else 0,
    }


def add_flags(rows):
    """Пометки: сравниваем каждый образец с медианой по всем и с порогами."""
    med_map = statistics.median(r["mapping_pct"] for r in rows)
    med_frag = statistics.median(r["frag_mean"] for r in rows)
    types = {r["lib_type"] for r in rows}

    # пол по парам: N3 и T3 → пациент "3"
    sex_by_patient = {}
    for r in rows:
        sex_by_patient.setdefault(r["sample"][1:], set()).add(r["sex"])

    for r in rows:
        f = []
        if r["mapping_pct"] < med_map - 10:      f.append("low_mapping")
        if r["frag_mean"] < 0.8 * med_frag:      f.append("short_fragments")
        if r["compat"] < 0.9:                    f.append("low_compat")
        if r["strand"] < 0.8:                    f.append("low_strandedness")
        if len(types) > 1:                       f.append("mixed_lib_types")
        if r["no_gene_tx"] > 0:                  f.append("tx_without_gene")
        if r["sex"] == "?":                      f.append("sex_unclear")
        if len(sex_by_patient[r["sample"][1:]] - {"?"}) > 1:
            f.append("sex_mismatch_in_pair")
        if r["MT_pct"] > 20:                     f.append("high_MT")
        r["flags"] = ",".join(f) if f else "ok"
    return med_map, med_frag


def main():
    if len(sys.argv) != 4:
        sys.exit("Использование: salmon_qc.py <папка salmon> <tx2gene.tsv> <выходной файл>")
    salmon_dir, tx2gene_path, out_path = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])

    tx2gene = load_tx2gene(tx2gene_path)
    sample_dirs = sorted((d for d in salmon_dir.iterdir() if (d / "quant.sf").exists()),
                         key=lambda d: patient_order(d.name))
    rows = [collect(d, tx2gene) for d in sample_dirs]
    med_map, med_frag = add_flags(rows)

    # формат чисел для таблицы
    fmt = {"processed_M": "{:.2f}", "mapped_M": "{:.2f}", "mapping_pct": "{:.1f}",
           "frag_mean": "{:.0f}", "compat": "{:.3f}", "strand": "{:.3f}",
           "XIST": "{:.0f}", "Y_genes": "{:.0f}", "ACTB": "{:.0f}", "GAPDH": "{:.0f}",
           "MT_pct": "{:.1f}", "top_gene_pct": "{:.1f}"}
    cols = list(rows[0].keys())

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        # шапка для MultiQC (custom content)
        f.write("# id: 'salmon_qc'\n")
        f.write("# section_name: 'Salmon: сводка проверок'\n")
        f.write(f"# description: 'Медиана mapping rate {med_map:.1f}%, длины фрагмента {med_frag:.0f} нт. "
                "Пары и фрагменты в млн; гены — сумма NumReads по tx2gene.'\n")
        f.write("# plot_type: 'table'\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(fmt.get(c, "{}").format(r[c]) for c in cols) + "\n")

    # короткая сводка в терминал
    print(f"Записано: {out_path}  ({len(rows)} образцов)")
    for r in rows:
        print(f"  {r['sample']:<5} mapping {r['mapping_pct']:5.1f}%  frag {r['frag_mean']:4.0f}  "
              f"strand {r['strand']:.2f}  sex {r['sex']}  → {r['flags']}")


if __name__ == "__main__":
    main()