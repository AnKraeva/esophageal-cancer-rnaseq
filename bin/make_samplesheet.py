#!/usr/bin/env python3
"""
Собирает samplesheet.csv для пайплайна по любому accession из ENA:
исследование (PRJNA..., SRP...), образец (SAMN...) или ран (SRR...).

Только стандартная библиотека Python — ничего устанавливать не нужно.

Сохраняет:
  metadata/ena_runs.tsv    — исходные метаданные ENA (запись о происхождении данных)
  metadata/samplesheet.csv — вход пайплайна: имя образца, ран, ссылки на FASTQ

"""
import argparse
import csv
import io
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
OUTDIR = Path("metadata")  # папка для метаданных; в отличие от data/, хранится в git
FIELDS = "run_accession,sample_alias,sample_title,library_layout,read_count,fastq_bytes,fastq_ftp"

parser = argparse.ArgumentParser(description="ENA accession -> samplesheet.csv")
parser.add_argument("accession", help="номер исследования, образца или рана")
args = parser.parse_args()

# --- Метаданные из ENA (ответ — TSV-текст) ----------------------------------
url = f"{ENA_API}?accession={args.accession}&result=read_run&fields={FIELDS}&format=tsv"
try:
    with urllib.request.urlopen(url, timeout=60) as resp:
        text = resp.read().decode()
except Exception as e:  # например, неверный номер или нет сети
    raise SystemExit(f"Не удалось получить данные для {args.accession}: {e}")

# Сохраняем ответ ENA как есть
OUTDIR.mkdir(exist_ok=True)
(OUTDIR / "ena_runs.tsv").write_text(text)

# Разбираем TSV в список словарей: [{"run_accession": "SRR...", ...}, ...]
runs = sorted(csv.DictReader(io.StringIO(text), delimiter="\t"),
              key=lambda r: r["run_accession"])
if not runs:
    raise SystemExit(f"Для {args.accession} в ENA не найдено ни одного рана")

# --- Samplesheet ------------------------------------------------------------
# Имя образца — sample_title из ENA (у нас "N21", "T4");
# если названия повторяются или пустые, берём номер рана — он уникален всегда
titles = [r["sample_title"] for r in runs]
use_titles = all(titles) and len(set(titles)) == len(titles)

with open(OUTDIR / "samplesheet.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample", "run_accession", "fastq_1", "fastq_2"])
    for r in runs:
        # В fastq_ftp ссылки через ";" (R1;R2); у single-end ссылка одна
        urls = ["https://" + u for u in r["fastq_ftp"].split(";")]
        fastq_2 = urls[1] if len(urls) > 1 else ""
        sample = r["sample_title"] if use_titles else r["run_accession"]
        writer.writerow([sample, r["run_accession"], urls[0], fastq_2])

layouts = sorted({r["library_layout"] for r in runs})
print(f"{OUTDIR}/samplesheet.csv: {len(runs)} образцов ({', '.join(layouts)})")