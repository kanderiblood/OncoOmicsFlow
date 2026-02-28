# Pipelines

## Paired Genomics
Содержит скрипт `exome_pipeline.sh` для обработки парных образцов (tumor/normal) в экзомных и WGS данных.

**Пример запуска:**
```bash
source configs/genomic_paired_config.sh
bash pipelines/paired_genomics/exome_pipeline.sh <TUMOR_SAMPLE> <NORMAL_SAMPLE>
Single Genomics (Onco Panel)

Содержит скрипт onco_pipeline.sh для одиночных образцов по таргетным онкопанелям.

Пример запуска:

source configs/onco_single_config.sh
bash pipelines/single_genomics/onco_pipeline.sh <SAMPLE>
RNA-seq

Содержит скрипты и конфиги для анализа RNA-секвенирования (обработка FastQ, выравнивание, QC, подсчет экспрессии и downstream-анализ).

Пример запуска:

source configs/rna_config.sh
bash pipelines/rna/rna_pipeline.sh <SAMPLE>
Proteomics

Содержит скрипты и конфиги для анализа протеомных данных (QuantMS и связанные инструменты).

Пример запуска:

source configs/proteomics_config.sh
bash pipelines/proteomics/proteomics_pipeline.sh
