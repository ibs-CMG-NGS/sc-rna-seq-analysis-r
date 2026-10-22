# [회신] expression_distribution pilot bundle 및 스키마 확인

- **수신**: cmg-seqviewer 개발팀
- **발신**: scRNA-seq 3차 분석 파이프라인 담당
- **관련 요청**: `EXPRESSION_DISTRIBUTION_ANALYSIS_TEAM_REQUEST.md`
- **상태**: pilot bundle 생성 완료 (per-sample pseudobulk)

---

## 1. 요약

- `dataset_type: "expression_distribution"` entry를 기존 `viewer_export` 번들에 추가했습니다.
- pilot은 **Dowsett10xnuc 1개 dataset**, **30개 marker genes**, **43개 cell type 전체**, **6개 biological sample**로 생성했습니다.
- 1순위 요청인 **per-sample pseudobulk**를 사용했습니다.
- 최종 `seqviewer_manifest.json`은 analysis 69 + distribution 1 = **70 datasets**입니다.
- 현재 viewer v1.4.0에는 `expression_distribution` 타입이 아직 없으므로,
  viewer 구현 후 pilot Import/UI 검증이 가능합니다.

## 2. pilot bundle

```text
output/tutorial_compare/export/seqviewer/
├── seqviewer_manifest.json
├── datasets/
│   ├── hypomap_Dowsett10xnuc_pseudobulk.parquet          # 7,740 rows, ~111 KB
│   ├── hypomap_Dowsett10xnuc_gene_index.parquet          # 30 genes
│   └── hypomap_Dowsett10xnuc_cell_type_summary.parquet   # 43 cell types
├── export_index.tsv
└── COLUMNS.md
```

- main table: `hypomap_Dowsett10xnuc_pseudobulk.parquet`
- row: `gene × cell_type × sample` (7,740 = 30 × 258)
- 258 sample × cell_type groups, 모두 `n_cells >= 10`
- `cell_type_summary`의 `n_cells` 합계 = 33,838 cells

## 3. sample / replicate 정의 (질문 #1)

pilot의 `sample`은 **biological sample/replicate**입니다.

- `dataset`: `Dowsett10xnuc`
- `sample`: `Dowsett10xnuc_adlib1`–`adlib4`, `Dowsett10xnuc_fasted1`–`fasted2`
- `condition`: `Normal chow` (4 samples), `Fasted` (2 samples)
- HypoMap 메타데이터에 별도 `donor` 컬럼은 없으며, 이 subset에서는 `sample_id`를 biological replicate로 사용했습니다.
- manifest에 `distribution_schema.sample_unit: "biological_replicate"`로 명시했습니다.

## 4. 확정 스키마

```json
{
  "schema_version": "1.0",
  "datasets": [
    {
      "dataset_id": "7ee0a902-4c61-8d94-4e1d-aa47c0246fb0",
      "alias": "HypoMap Dowsett10xnuc pseudobulk distribution",
      "file_path": "hypomap_Dowsett10xnuc_pseudobulk.parquet",
      "dataset_type": "expression_distribution",
      "organism": "Mus musculus",
      "row_count": 7740,
      "gene_count": 30,
      "cell_type_count": 43,
      "sample_count": 6,
      "distribution_schema": {
        "schema_version": "1.0",
        "layout": "pseudobulk",
        "gene_id_type": "symbol",
        "cell_type_column": "cell_type",
        "sample_column": "sample",
        "sample_unit": "biological_replicate",
        "dataset_column": "dataset",
        "condition_column": "condition",
        "group_columns": ["condition", "dataset"],
        "value_column": "mean_expression",
        "pct_column": "pct_expressing",
        "n_cells_column": "n_cells",
        "pct_scale": "0-100",
        "detection": "raw counts > 0",
        "value_definition": "mean of normalized value layer per sample x cell_type",
        "n_genes": 30,
        "n_samples": 6,
        "min_cells": 10,
        "companion_files": {
          "gene_index": "hypomap_Dowsett10xnuc_gene_index.parquet",
          "cell_type_summary": "hypomap_Dowsett10xnuc_cell_type_summary.parquet"
        }
      }
    }
  ]
}
```

parquet columns:

```text
gene_id, gene_name, cell_type, sample, dataset, condition,
n_cells, n_expressing, pct_expressing, mean_expression, mean_counts
```

- `mean_expression`: normalized `data` layer의 sample × cell_type 평균
- `mean_counts`: raw counts 평균
- `pct_expressing`: `100 * n_expressing / n_cells`
- cell-level 값은 포함하지 않습니다.

## 5. 요청 질문 답변

| # | 질문 | 답변 |
|---|---|---|
| 1 | HypoMap sample/replicate 정의 | Dowsett10xnuc subset은 6 biological samples(4 Normal chow, 2 Fasted)입니다. `Dataset`은 study, `sample_id`가 replicate입니다. |
| 2 | per-sample pseudobulk export 가능? | 가능합니다. pilot이 per-sample pseudobulk입니다. |
| 3 | quantile summary 가능? | 가능하지만 pilot에는 포함하지 않았습니다. 필요 시 p10/p25/median/p75/p90 fallback을 추가할 수 있습니다. |
| 4 | 선택 gene cell-level export 가능? | 가능합니다. pilot에는 포함하지 않았고, gene-partitioned parquet 또는 long parquet로 별도 협의 후 제공할 수 있습니다. |
| 5 | gene 범위 | pilot은 30개 marker genes입니다. 전체 57k는 viewer의 lazy query 전략과 전달 크기를 협의한 뒤 결정하는 것이 좋습니다. |
| 6 | condition/sample/dataset/donor 메타데이터 | `sample`, `dataset`, `condition` 포함. HypoMap에는 별도 `donor` 컬럼이 없습니다. |
| 7 | pilot bundle 제공 가능? | 제공 완료했습니다. 위 경로 참고. |
| 8 | `viewer_export` manifest merge 제약? | 없습니다. merge 스크립트에 distribution staging을 추가해 최종 manifest에 1 entry로 병합합니다. |
| 9 | 조건 비교/통계 목적? | pilot은 descriptive violin/box용입니다. 조건 간 통계 검정이 목적이면 pseudobulk DE 결과를 별도로 사용해야 합니다. |
| 10 | pilot 파일 크기/전달 위치 | main parquet ~111 KB, companions 포함 ~124 KB. 로컬 `output/tutorial_compare/export/seqviewer/`에 있으며 필요 시 업로드/전달하겠습니다. |

## 6. 검증

- `Agrp` / `C66-46: Agrp.GABA-4` / `Dowsett10xnuc_fasted1`에 대해 원본 Seurat에서 독립 재계산:
  - `n_cells` 116 = 116
  - `mean_expression` 1.993835 = 1.993835
  - `pct_expressing` 91.37931 = 91.37931
  - `mean_counts` 6.525862 = 6.525862
- `cell_type_summary`의 `n_cells` 합계 = 33,838 (subset 전체 세포 수)
- 최종 manifest: 70 datasets (analysis 69 + distribution 1), companion 파일 존재 확인

## 7. 다음 단계

1. viewer에서 pilot bundle의 violin/box UI·성능 검증
2. `distribution_schema` 필드 최종 확정
3. 전체 atlas 범위, gene panel, quantile/cell-level 옵션 협의
