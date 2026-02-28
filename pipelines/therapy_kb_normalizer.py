# therapy_kb_normalizer.py
"""
Модуль нормализации и унификации записей из разных баз знаний терапии.
Цель: привести данные к единой структуре, независимо от источника.

Поддерживаемые источники (на старте):
- OncoKB
- CIViC
- COSMIC (Cancer Gene Census + actionable variants)
- (можно легко добавить JAX CKB, MolecularMatch, etc.)

Единая структура записи (TherapyRecord):
{
    "source": str,                  # "oncokb", "civic", "cosmic", ...
    "gene": str,                    # HGNC symbol, нормализованный
    "variant": str,                 # HGVS на уровне белка или cDNA, предпочтительно protein
    "variant_type": str,            # SNV, indel, fusion, CNV, etc.
    "tumor_type": list[str],        # онтология или свободный текст
    "therapy": list[str],           # названия препаратов / комбинаций
    "evidence_level": str,          # A, B, C, D, R1, etc. или числовой эквивалент
    "sensitivity": str,             # "sensitive", "resistant", "intermediate", ...
    "direction": str,               # "sensitive", "resistant", "no evidence", ...
    "pmids": list[int],             # PMID источников
    "description": str,             # краткое описание утверждения
    "last_updated": str or None,    # дата обновления записи
    "original_record": dict,        # сырые данные источника (для отладки)
}
"""

from typing import Dict, List, Any, Optional
import re
from datetime import datetime


class TherapyRecord:
    """Единая нормализованная запись о терапевтической значимости варианта"""

    def __init__(self):
        self.source: str = ""
        self.gene: str = ""
        self.variant: str = ""
        self.variant_type: str = ""
        self.tumor_type: List[str] = []
        self.therapy: List[str] = []
        self.evidence_level: str = ""
        self.sensitivity: str = ""
        self.direction: str = ""
        self.pmids: List[int] = []
        self.description: str = ""
        self.last_updated: Optional[str] = None
        self.original_record: Dict[str, Any] = {}

    def to_dict(self) -> Dict[str, Any]:
        return {
            k: v for k, v in self.__dict__.items()
            if not k.startswith("_") and v is not None and v != []
        }

    def __repr__(self):
        return f"TherapyRecord(gene={self.gene}, variant={self.variant}, source={self.source})"


class TherapyNormalizer:
    """Основной класс нормализации"""

    def __init__(self):
        self.source_handlers = {
            "oncokb": self._normalize_oncokb,
            "civic": self._normalize_civic,
            "cosmic": self._normalize_cosmic,
            # Добавляйте новые источники сюда
        }

    def normalize(self, record: Dict[str, Any], source: str) -> Optional[TherapyRecord]:
        """
        Главный метод: принимает сырую запись и имя источника → возвращает нормализованную запись
        """
        source = source.lower().strip()
        handler = self.source_handlers.get(source)
        if not handler:
            print(f"Неизвестный источник: {source}")
            return None

        try:
            normalized = handler(record)
            if normalized:
                normalized.source = source
                normalized.original_record = record
                return normalized
            else:
                return None
        except Exception as e:
            print(f"Ошибка нормализации записи из {source}: {e}")
            return None

    # oncokb
    def _normalize_oncokb(self, record: Dict) -> Optional[TherapyRecord]:
        tr = TherapyRecord()

        tr.gene = record.get("hugoSymbol", "").strip().upper()
        tr.variant = record.get("alteration", "").strip()

        # Вариант типа
        tr.variant_type = record.get("variantType", "SNV").strip()

        # Опухоль
        if "cancerTypes" in record:
            tr.tumor_type = [ct["name"] for ct in record["cancerTypes"] if ct.get("name")]

        # Терапия
        if "drugs" in record:
            tr.therapy = [d["drugName"] for d in record["drugs"] if d.get("drugName")]

        tr.evidence_level = record.get("level", "")
        tr.sensitivity = record.get("oncogenic", "").lower()
        tr.direction = "sensitive" if "sensitive" in tr.sensitivity else "resistant"

        # PMID
        if "pmids" in record:
            tr.pmids = [int(p) for p in record["pmids"] if p.isdigit()]

        tr.description = record.get("abstract", "") or record.get("description", "")

        if "lastUpdate" in record:
            tr.last_updated = record["lastUpdate"]

        return tr if tr.gene and tr.variant else None
    # civic 
    def _normalize_civic(self, record: Dict) -> Optional[TherapyRecord]:
        tr = TherapyRecord()

        tr.gene = record.get("gene", {}).get("name", "").strip().upper()

        variant_data = record.get("variant", {})
        tr.variant = variant_data.get("name", "").strip()
        tr.variant_type = variant_data.get("variant_types", [{}])[0].get("name", "unknown")

        # Опухоль / заболевание
        disease = record.get("disease", {})
        if disease and disease.get("name"):
            tr.tumor_type = [disease["name"]]

        # Терапия
        if "therapies" in record:
            tr.therapy = [t.get("name", "") for t in record["therapies"] if t.get("name")]

        # Уровень доказательности
        evidence = record.get("evidence_items", [{}])[0] if record.get("evidence_items") else {}
        tr.evidence_level = evidence.get("evidence_level", "")
        tr.direction = evidence.get("clinical_significance", "").lower()
        tr.sensitivity = tr.direction  # часто совпадает

        if "pubmed_id" in evidence:
            try:
                tr.pmids = [int(evidence["pubmed_id"])]
            except:
                pass

        tr.description = evidence.get("description", "")

        if "updated_at" in record:
            tr.last_updated = record["updated_at"]

        return tr if tr.gene and tr.variant else None

    #cosmic
    def _normalize_cosmic(self, record: Dict) -> Optional[TherapyRecord]:
        tr = TherapyRecord()

        tr.gene = record.get("gene_symbol", "").strip().upper()
        tr.variant = record.get("mutation_aa", "") or record.get("cds_mutation", "")

        tr.variant_type = record.get("mutation_type", "SNV")

        tr.tumor_type = [record.get("cancer_type", "")] if record.get("cancer_type") else []

        # Терапия — в COSMIC часто косвенно
        tr.therapy = [record.get("drug", "")] if record.get("drug") else []
        tr.evidence_level = record.get("tier", "") or "COSMIC"
        tr.direction = record.get("effect", "").lower()

        if "pubmed" in record:
            try:
                tr.pmids = [int(record["pubmed"])]
            except:
                pass

        tr.description = record.get("description", "")

        return tr if tr.gene and tr.variant else None


#пример для использования 

if __name__ == "__main__":
    normalizer = TherapyNormalizer()

    # Пример сырой записи OncoKB
    oncokb_example = {
        "hugoSymbol": "BRAF",
        "alteration": "V600E",
        "variantType": "missense",
        "cancerTypes": [{"name": "Melanoma"}],
        "drugs": [{"drugName": "Vemurafenib"}],
        "level": "A",
        "pmids": ["12345678", "87654321"],
        "lastUpdate": "2025-11-15"
    }

    result = normalizer.normalize(oncokb_example, "oncokb")
    if result:
        print("Нормализованная запись:")
        import json
        print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))
