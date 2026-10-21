# cmg-seqviewer 연동 문의 (scRNA-seq 3차 분석 파이프라인 → 뷰어)

> 발신: scRNA-seq 3차 분석 파이프라인(R/Seurat) 담당
> 수신: cmg-seqviewer 개발팀
> 목적: 본 파이프라인 산출물을 cmg-seqviewer로 넘겨 연구자가 검색·비교·시각화할 수 있게 하는 연동 스펙 확정

---

## 1. 배경

- 본 파이프라인: `ibs-CMG-NGS/sc-rna-seq-analysis-r` (R/Seurat + Snakemake, **3차 분석 전용**)
  — 기존 `sc-rna-seq-analysis`(Scanpy판)과는 **별개 구현**입니다.
- 연동 원칙: **뷰어는 Seurat/h5ad를 직접 읽지 않습니다.** 본 파이프라인이 뷰어의 기존 스키마(DE / GO·KEGG / Multi-Group)로 reshape한 **parquet + manifest**를 `output/<run>/export/seqviewer/`에 내보내고, 뷰어가 이를 소비합니다.
- 기대 효과: 세포타입별 DE/GO 결과를 DB에 적재 → **Cell Type / Organism / Tags 검색** + **Cross-Dataset Comparison**(세포타입 간 비교) + Volcano/Dot plot 등 기존 기능 재사용.

## 2. 확인한 사항 (문서 기준)

- DE: `gene_id`, `log2fc`, `adj_pvalue`, `baseMean` (+pvalue)
- GO/KEGG: `term_id`, `description`, `gene_count`, `fdr`, `gene_symbols`, `gene_set`, `ontology`, `gene_ratio`
- ATAC: `peak_id` 기반 자동인식
- CLI: `cmg-seqviewer-cli load/filter/run/batch` (headless)
- DB: `database/datasets/*.parquet` + `metadata.json`, v1.2.4에서 `seqviewer_manifest.json` 인식 추가

---

## 3. 문의 (우선순위순)

### Q1. [높음] `seqviewer_manifest.json` 정식 스키마
- **현황**: v1.2.4 릴리즈노트에 "`seqviewer_manifest.json` recognized during database import" 언급만 있고, 스키마 문서가 없습니다.
- **필요 정보**:
  1. 허용 필드 목록 · 각 필드 필수 여부 · 값 형식 · 최소 예시 JSON
  2. DB 자동인식 폴더명 규칙 `/YYYY-initials-organism-type/`의 정확한 파싱 규칙(manifest 부재 시 추론 로직 포함)
- **이유**: 본 파이프라인이 manifest를 생성해 DB에 적재 → 연구자가 세포타입/조건/태그로 검색하는 흐름을 구현하려 합니다.
- **우리 계획 필드**: `alias`, `dataset_type`, `cell_type`, `organism`, `tissue`, `condition`, `timepoint`, `tags`, `notes`, `row_count`, `gene_count`, `significant_genes`, `provenance`(파이프라인·annotation 출처·버전), `schema_version`

### Q2. [높음] Multi-Group dataset의 정확한 요구 스키마
- **현황**: 릴리즈노트에 "Multi-Group dataset detection recognizes `cluster_id` + numeric sample columns, no `padj`" 로만 기술되어 있습니다.
- **필요 정보**:
  1. gene id 첫 컬럼의 표준명(자동인식되는 헤더 목록)
  2. numeric sample 컬럼의 명명/타입 규칙(예: `S1, S2, …`)
  3. `cluster_id`가 **없어도** Multi-Group으로 인식되는지(단순 발현행렬로 PCA/Heatmap 가능 여부) vs 필수인지
  4. DE로 오탐되지 않고 Multi-Group으로 판정되는 **정확한 조건**
- **이유**: pseudobulk count matrix(gene × `cell_type__sample`)를 그대로 넘겨 **세포타입별 PCA / Heatmap / replicate bar+scatter**를 쓰려 합니다.

### Q3. [중] DE dataset에 추가 메타 컬럼 유지 가능 여부
- 본 파이프라인 DE 출력에는 `cell_type`, `contrast` 컬럼이 이미 포함됩니다. 이들을 **그대로 남겨서** 넘겨도 자동인식·필터·시각화에 충돌이 없습니까? (제거 후 넘기는 것이 권장이면 알려주세요.)

### Q4. [중] CLI `load`의 기계 판독 출력
- `cmg-seqviewer-cli load <file>`이 판정 결과(dataset_type, rows, 필수 컬럼 충족 여부)를 **JSON 등으로 출력하는 옵션**이 있습니까? exit code 규약은?
- **이유**: 본 파이프라인의 CI에서 export 파일이 올바른 type으로 인식되는지 **자동 검증**하고 싶습니다.

### Q5. [중] GO/KEGG `Gene Set` / `Ontology` 명시 컬럼
- 세포타입×contrast×방향×ontology 조합이 많아 **시트명 추론으로는 부족**합니다. `Gene Set`/`Ontology`를 **명시 컬럼**으로 주는 방식을 권장하십니까? (문서상 `✔*` 로 보이는데 확정 부탁드립니다.)

### Q6. [낮음] 다수 dataset의 DB browser 성능
- 세포타입×contrast 조합으로 dataset 수가 수십~수백 개가 될 수 있습니다. DB browser의 **그룹화/폴더 로드** 기능 일정이 있습니까?

### Q7. [낮음] roadmap v2.0 single-cell 지원
- v2.0에 "Single-cell RNA-Seq support"가 예정되어 있습니다. 계획 중인 스키마를 공유해 주시면 본 파이프라인 export를 그 방향에 미리 정렬하겠습니다.

---

## 4. 기능 요청 (선택)

### R1. `composition` dataset type
- 본 파이프라인은 **cell type 비율 변화(propeller, limma 기반)** 결과를 산출합니다(`cell_type`, `PropMean.<cond>`, `PropRatio`, `P.Value`, `FDR`).
- 현재 뷰어에 대응하는 dataset type / plot이 없습니다. **신규 `composition` type + 전용 plot**을 제안합니다.
- 우회안: cell_type을 gene처럼, 조건별 비율을 sample 컬럼처럼 넣어 **Multi-Group heatmap**으로 표현 가능한지도 확인 부탁드립니다.

---

## 5. 다음 단계 (회신 후)

| 단계 | 내용 | 의존 |
| --- | --- | --- |
| E1 | DE / GO·KEGG 어댑터 export (뷰어 무수정) | — (독립 진행) |
| E2 | `manifest.json` 생성 + DB 적재 | Q1 |
| E3 | Multi-Group export + upstream `cluster_id` | Q2 |
| E4 | composition type 등 신규 기능 반영 | R1 |

회신 주시면 E2·E3 설계에 즉시 반영하겠습니다. 감사합니다.
