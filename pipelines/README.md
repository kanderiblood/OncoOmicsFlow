# Pipelines

## Paired Genomics
Содержит скрипт `exome_pipeline.sh` для обработки парных образцов (tumor/normal) в WGS,WES,target

**Пример запуска:**
```bash
source configs/genomic_paired_config.sh
bash pipelines/paired_genomics/exome_pipeline.sh <TUMOR_SAMPLE> <NORMAL_SAMPLE>
##Single Genomics 

Содержит скрипт onco_pipeline.sh для одиночных образцов 

Пример запуска:

source configs/onco_single_config.sh
bash pipelines/single_genomics/onco_pipeline.sh <SAMPLE>
##RNA-seq

Содержит скрипты для анализа RNA-секвенирования 

Пример запуска:

source configs/rna_config.sh
bash pipelines/rna/rna_pipeline.sh <SAMPLE>
##Proteomics

Содержит скрипты анализа протеомных данных 

Пример запуска:

source configs/proteomics_config.sh
bash pipelines/proteomics/proteomics_pipeline.sh
