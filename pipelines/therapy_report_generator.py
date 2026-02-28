# therapy_report_generator.py
"""
Модуль автоматического формирования клинических отчётов на основе унифицированной базы знаний терапии.
"""

from typing import List, Dict, Any, Optional
from therapy_kb_normalizer import TherapyRecord, TherapyNormalizer   # ← импортируем из предыдущего файла
import datetime


class TherapyKnowledgeBase:
    """Хранилище нормализованных записей (in-memory база)"""

    def __init__(self):
        self.records: List[TherapyRecord] = []

    def add_record(self, record: TherapyRecord):
        self.records.append(record)

    def add_records(self, records: List[TherapyRecord]):
        self.records.extend(records)

    def query(self,
              gene: Optional[str] = None,
              variant: Optional[str] = None,
              tumor_type: Optional[str] = None,
              min_evidence: Optional[str] = None) -> List[TherapyRecord]:
        """Простой поиск по базе"""
        result = []
        for r in self.records:
            if gene and r.gene.upper() != gene.upper():
                continue
            if variant and variant.upper() not in r.variant.upper():
                continue
            if tumor_type and not any(tumor_type.lower() in t.lower() for t in r.tumor_type):
                continue

            # фильтр по уровню доказательности (A > B > C > D)
            if min_evidence and self._evidence_score(r.evidence_level) < self._evidence_score(min_evidence):
                continue

            result.append(r)
        return sorted(result, key=lambda x: self._evidence_score(x.evidence_level), reverse=True)

    @staticmethod
    def _evidence_score(level: str) -> int:
        mapping = {"A": 4, "B": 3, "C": 2, "D": 1, "R1": 3, "R2": 2, "": 0}
        return mapping.get(level.upper(), 0)


class ClinicalReportGenerator:
    """Генератор персонализированных клинических отчётов"""

    def __init__(self, kb: TherapyKnowledgeBase):
        self.kb = kb

    def generate_markdown_report(self,
                                 patient_id: str,
                                 tumor_type: str,
                                 patient_variants: List[Dict[str, str]],   # [{"gene": "BRAF", "variant": "V600E"}, ...]
                                 title: str = "Персонализированный клинико-геномный отчёт") -> str:

        report = []
        report.append(f"# {title}")
        report.append(f"**Пациент:** {patient_id}")
        report.append(f"**Тип опухоли:** {tumor_type}")
        report.append(f"**Дата отчёта:** {datetime.date.today().isoformat()}\n")

        report.append("## Выявленные варианты")
        for v in patient_variants:
            report.append(f"- **{v['gene']} {v['variant']}**")

        report.append("\n## Рекомендуемые терапевтические опции\n")

        found_any = False
        for var in patient_variants:
            matches = self.kb.query(
                gene=var["gene"],
                variant=var["variant"],
                tumor_type=tumor_type,
                min_evidence="C"          # минимальный уровень C и выше
            )

            if matches:
                found_any = True
                report.append(f"### {var['gene']} {var['variant']}")
                for m in matches[:8]:     # топ-8 самых сильных доказательств
                    therapies = ", ".join(m.therapy) if m.therapy else "—"
                    report.append(
                        f"- **{therapies}**  \n"
                        f"  Уровень доказательности: **{m.evidence_level}**  \n"
                        f"  Эффект: {m.direction.capitalize()}  \n"
                        f"  Тип опухоли: {', '.join(m.tumor_type)}  \n"
                        f"  Источники: {m.source.upper()}"
                    )
                report.append("")

        if not found_any:
            report.append("**Не найдено клинически значимых совпадений в базе знаний.**\n")

        report.append("## Примечания")
        report.append("- Отчёт сформирован автоматически на основе унифицированной базы OncoOmicsFlow.")
        report.append("- Рекомендации носят информационный характер и требуют подтверждения клиническим специалистом.")
        report.append(f"- Всего записей в базе: {len(self.kb.records)}")

        return "\n".join(report)

#пример использования 

if __name__ == "__main__":
    # Создаём базу
    kb = TherapyKnowledgeBase()
    normalizer = TherapyNormalizer()

    # Пример добавления записей (можно загружать из JSON/CSV)
    examples = [
        {"hugoSymbol": "BRAF", "alteration": "V600E", "cancerTypes": [{"name": "Melanoma"}], "drugs": [{"drugName": "Vemurafenib"}], "level": "A", "source": "oncokb"},
        {"gene": {"name": "EGFR"}, "variant": {"name": "L858R"}, "disease": {"name": "Non-small cell lung cancer"}, "therapies": [{"name": "Osimertinib"}], "evidence_level": "A", "source": "civic"},
    ]

    for ex in examples:
        src = ex.pop("source")
        rec = normalizer.normalize(ex, src)
        if rec:
            kb.add_record(rec)

    # Генерируем отчёт
    generator = ClinicalReportGenerator(kb)

    patient_data = [
        {"gene": "BRAF", "variant": "V600E"},
        {"gene": "EGFR", "variant": "L858R"}
    ]

    markdown_report = generator.generate_markdown_report(
        patient_id="PAT-2025-001",
        tumor_type="Melanoma",
        patient_variants=patient_data
    )

    print(markdown_report)

    # Сохраняем в файл
    with open("clinical_report_PAT-2025-001.md", "w", encoding="utf-8") as f:
        f.write(markdown_report)

    print("\nОтчёт сохранён в clinical_report_PAT-2025-001.md")
