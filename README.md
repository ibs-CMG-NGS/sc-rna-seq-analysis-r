# scRNA-seq Downstream Analysis Pipeline

R과 Snakemake 기반의 재현 가능한 단일세포 RNA-seq 3차 분석 파이프라인입니다.  
Seurat (`.rds`) 또는 AnnData (`.h5ad`) 형식을 입력으로 받습니다. 목적에 따라 **(1) 반복 없이 유전자 발현을 요약하는 탐색 분석**과 **(2) 생물학적 반복을 이용하는 비교 통계 분석**으로 실행할 수 있으며, 비교 분석은 세포 조성 → Pseudobulk DE → 유전자 세트 풍부화 → (선택) Trajectory → (선택) 세포 간 커뮤니케이션 → 통합 HTML 보고서까지 자동 생성합니다.

---

## 주요 기능

| Phase | 모듈 | 주요 내용 |
|:---:|---|---|
| 1 | **입력 처리** | Seurat / AnnData 자동 변환, 메타데이터 검증, 그룹 크기 필터 |
| 2a | **Pseudobulk DE** | DESeq2 또는 edgeR, 공변량 지원, Volcano plot 자동 생성 |
| 2b | **Compositional 분석** | propeller (speckle)로 세포 비율 차이 검정 |
| 3 | **Enrichment** | GO (BP/CC/MF) + KEGG, clusterProfiler ORA, dotplot 자동 생성 |
| 4 | **Trajectory** *(선택)* | Monocle3 또는 Slingshot, 생물학적 판단에 기반한 root cell 지정 |
| 5 | **Cell-Cell Communication** *(선택)* | CellChat 또는 LIANA, 조건 간 비교 지원 |
| 6 | **보고서** | R Markdown HTML 보고서, Methods 초안 (`methods_draft.md`) 자동 생성 |
| 7 | **cmg-seqviewer 연동** | DE / GO·KEGG / Composition / Multi-Group 결과를 뷰어 Import Folder 형식(parquet + `seqviewer_manifest.json`)으로 export |
| 8 | **Canonical tables** | 분석·발현 결과를 `tables/*.csv` + `export/excel/*.xlsx`로도 export |

## 분석 모드: 탐색 vs 비교(통계)

이 파이프라인은 목적이 다른 두 분석을 같은 config 구조로 실행합니다. **가장 흔한 실수는 생물학적 반복이 없는 데이터에 비교 통계 모듈을 켜는 것입니다.**

| | 탐색·서술 분석 | 비교·통계 분석 |
|---|---|---|
| 질문 | "이 유전자가 어느 세포타입에서 발현되나?" | "조건 간 발현·세포 조성·경로가 다른가?" |
| 필요 데이터 | annotated Seurat / AnnData | + `sample_id`(생물학적 샘플)와 `condition` |
| 최소 반복 | 없음 | 조건당 생물학적 샘플 **≥2** (권장 ≥3) |
| 핵심 config | `optional_modules.exploration: true`, `modules.*: false` | `modules.pseudobulk_de/composition/enrichment: true`, `optional_modules.exploration: false` |
| 주요 출력 | `exploration/*.csv`, bar/dot/heatmap | `pseudobulk_de/*`, `composition/*`, `enrichment/*`, volcano plot, `tables/`, `export/excel/` |
| 통계 | 없음 (서술적) | DESeq2/edgeR, propeller, GO/KEGG ORA |
| 주의 | `split_by: Dataset`으로 dataset 교란 확인 | 세포를 반복으로 세면 pseudoreplication |

- **탐색 분석**은 `exploration` 모듈만 켭니다. 실제 비교 설계가 없어도 되지만, 현재 loader는 `cell_type`/`sample`/`condition` 컬럼 자체를 요구하므로 더미 값(예: `condition` = `"all"`)을 지정해야 합니다. p-value는 계산하지 않습니다.
- **비교 분석**은 `sample_column`이 **세포 ID가 아닌 생물학적 샘플 ID**여야 합니다. 조건당 샘플이 1개면 재현성 있는 통계 검정이 어렵습니다.
- **둘 다** 켤 수 있습니다. 이 경우 탐색 결과와 비교 결과가 같은 HTML 리포트에 포함됩니다.
- `modules.enrichment: true`는 `modules.pseudobulk_de: true`가 필요합니다.
- 바로 복사해 쓸 수 있는 시작 파일: [config_explore.yml](configs/template/config_explore.yml), [config_compare.yml](configs/template/config_compare.yml)
- 실행 예시와 실제 출력은 [HypoMap 튜토리얼](docs/tutorial_hypomap.md)의 **5장(탐색)**과 **6장(비교)**을 참고하세요.

### 탐색 전용 최소 설정

```yaml
single_cell:
  input_path: "path/to/annotated.rds"
  input_format: "seurat"
  cell_type_column: "cell_type"
  sample_column: "sample_id"
  condition_column: "condition"
  cell_types: null

modules:
  pseudobulk_de: false
  composition: false
  enrichment: false

optional_modules:
  exploration: true

exploration:
  enabled: true
  genes: ["Agrp", "Pomc", "Glp1r"]
  group_by: "cell_type"
  split_by: "Dataset"
  min_cells: 10
```

### 비교 통계 최소 설정

```yaml
single_cell:
  input_path: "path/to/annotated.rds"
  input_format: "seurat"
  cell_type_column: "cell_type"
  sample_column: "sample_id"      # 세포 ID가 아닌 생물학적 샘플 ID
  condition_column: "condition"
  cell_types: null

modules:
  pseudobulk_de: true
  composition: true
  enrichment: true

de_analysis:
  method: "DESeq2"
  reference_level: "control"
```

## cmg-seqviewer 연동 (export)

`viewer_export` 하나로 분석 결과와 탐색용 발현 요약을 같은 Import Folder에 내보냅니다.

```yaml
viewer_export:
  enabled: true
  analysis: true            # DE / GO·KEGG / composition / Multi-Group
  expression: true          # 전체 gene x cell_type 평균 발현 요약
  schema_version: "1.0"
  expression_prefix: "my_experiment"
  expression_alias: "My experiment mean expression"
  expression_cell_type_level: "C185"
  expression_primary: "wide"        # wide 권장; long은 per-dataset (viewer v1 5M rows 한도)
  expression_include_long: false    # wide일 때 long 테이블도 추가
  expression_pct_scale: "0-100"
```

- `analysis: true`는 `modules.pseudobulk_de`가 켜져 있어야 합니다.
- `expression: true`는 `expression_export` 설정을 사용해 전체 발현 요약을 자동 실행합니다.
- 두 결과는 내부에서 별도로 생성된 뒤 **하나의 manifest**로 병합됩니다.
- 출력은 항상 `output/<실험명>/export/seqviewer/` 한 곳입니다.

```text
export/seqviewer/
├── seqviewer_manifest.json      # 모든 데이터셋 목록
├── datasets/*.parquet
├── export_index.tsv             # 파일별 유형·세포타입·조건 목록
├── COLUMNS.md                   # 컬럼 규약
└── export_complete.flag
```

| 파일 | dataset_type | 내용 |
|---|---|---|
| `DE__<cell_type>__<contrast>.parquet` | differential_expression | 세포타입 × contrast별 DE 결과 |
| `ENR__<cell_type>_<contrast>__GO_<ONT>_<up\|down>.parquet`, KEGG 동일 패턴 | go_analysis | GO/KEGG 결과 (세포타입 × contrast × 방향) |
| `COMP__propeller.parquet` | composition | 세포타입 비율 검정 (propeller) |
| `MG__<cell_type>.parquet` | multi_group | pseudobulk raw count, 샘플 열은 `<condition>_<rep>` |
| `<prefix>_wide_mean.parquet` (wide) 또는 `<prefix>_long.parquet` (long) | expression_summary | 전체 gene × cell_type 평균 발현 |
| `<prefix>_wide_pct.parquet`, `<prefix>_gene_index.parquet`, `<prefix>_cell_type_summary.parquet` | companion | layout에 따라 일부만 생성되는 부속 파일 |

뷰어에서는 `output/<실험명>/export/seqviewer/` 폴더 자체를 **Import Folder**로 선택합니다. `metadata.json`이 아니라 `seqviewer_manifest.json`을 읽습니다.

> **viewer 호환성 (2026-10-21 확인)**: cmg-seqviewer `feat/cli` `46a491a`(v1.4.0 예정)에서 `composition`과 `expression_summary`를 지원합니다. `tutorial_compare`는 69/69, C185/C66/long 발현 번들은 각 1개씩 정상 Import되는 것을 viewer팀이 실측했습니다. v1.3.2 이하에서는 두 타입이 skip되므로 새 번들을 보려면 viewer 업그레이드가 필요합니다.
>
> C66 legacy `multi_group` fallback은 v1.4.0 배포 후 제거 예정입니다. fallback과 `expression_summary`를 같은 DB에 동시에 Import하면 중복 dataset으로 보이므로 하나만 사용하세요.
>
> 기존 `export_seqviewer` / `expression_export.enabled` 설정은 하위 호환을 위해 계속 동작합니다. 새 config는 `viewer_export` 사용을 권장합니다.

## Canonical CSV / Excel tables

`viewer_export`와 별개로, 분석·발현 결과를 **CSV 원본 테이블**로도 남깁니다.
viewer parquet은 서비스용 어댑터이고, `tables/` 아래 CSV가 재사용 가능한
canonical 결과입니다.

```yaml
table_export:
  enabled: true       # output/<실험명>/tables/ 생성+병합
  excel: true         # output/<실험명>/export/excel/*.xlsx 생성
```

- `tables/manifest.json`에 모든 CSV의 경로·유형·세포타입·조건이 기록됩니다.
- `tables/index.tsv`로 파일 목록을 빠르게 확인할 수 있습니다.
- Excel은 `dataset_type`별 workbook으로 묶입니다:
  `differential_expression.xlsx`, `go_analysis.xlsx`, `composition.xlsx`,
  `multi_group.xlsx`, `expression_summary.xlsx`.
- Excel 한계(1,048,575 data rows / 16,383 columns)를 넘는 표는 CSV-only로
  남고 `export/excel/excel_manifest.json`의 `skipped`에 기록됩니다.

```text
tables/
├── manifest.json
├── index.tsv
├── README.md
├── analysis/datasets/*.csv
└── expression/datasets/*.csv
export/excel/
├── differential_expression.xlsx
├── go_analysis.xlsx
├── composition.xlsx
├── multi_group.xlsx
├── expression_summary.xlsx
├── excel_manifest.json
└── README.md
```

---

### Phase 0a — Label Transfer (optional, reference-map)

`cell_type` 라벨이 없는 외부 데이터(예: HypoMap 같은 published scRNA-seq)를 위해 frozen reference model에 매핑하여 라벨을 보조한다 (`mapscvi` 패키지 사용). 설정: `optional_modules.label_transfer: true` 시 `load_input` 이전에 실행되며, 소요 시간은 세포 수·GPU 여부에 따라 증가한다.

| 키 | 설명 |
|---|---|
| `label_transfer.method` | 현재 `mapscvi`만 지원 (HypoMap frozen scVI model) |
| `label_transfer.label_column` | reference meta column (`C66_named` 또는 `C185_named`) |
| `label_transfer.max_epochs` | scVI mapping epochs (10–20 권장) |

<details><summary>설치·실행 안내 (클릭)</summary>

```bash
# 1) conda 환경에 scvi-tools 포함 후 활성화 (environment.yml 참조)
conda activate scrnaseq-analysis
# 2) mapscvi 설치 (R에서 remotes를 사용)
Rscript scripts/install_mapscvi.R
# 3) config에서 optional_modules.label_transfer: true 설정 후 실행
```

출력: `output/<실험명>/label_transfer/` 아래 `annotated.rds` (predicted + prediction_probability 반영, `cell_type` 컬럼이 없었던 경우에만 predictions로 채움) + `label_provenance.csv` (세포별 라벨·확률) + `label_provenance.yml` (모델 경로, 패키지 버전, input 형태 등 재현성 정보).
</details>

**bulk RNA-seq 파이프라인과의 일관성 유지:**
- `configs/template/config.yml` (Git 추적) + `configs/config_*.yml` (gitignore)
- `padj_cutoff`, `log2fc_cutoff` 최상위 레벨 키
- `profiles/default` Snakemake profile 기반 direct `snakemake` 실행
- `run_pipeline.sh` / `run_downstream_only.sh`는 deprecated 호환 shim
- `config_used.yml` 자동 저장, 출력물 CSV + Excel + PNG

---

## 시작하기

> **처음 사용하는 분께**: 환경 설치부터 HypoMap 다운로드, 탐색적 분석, 비교 통계 분석까지
> 단계별로 설명한 [HypoMap 튜토리얼](docs/tutorial_hypomap.md)을 먼저 따라 해 보세요.

### 1. 저장소 클론

```bash
git clone <repo_url>
cd scran-seq-analysis
```

### 2. Conda 환경 설치

```bash
conda env create -f environment.yml
conda activate scrnaseq-analysis
```

> **채널:** `conda-forge` → `bioconda` 순서. `defaults` 채널 미포함 (Anaconda ToS 자유).

### 3. Config 파일 작성

먼저 분석 모드에 맞는 **시작 템플릿**을 복사합니다.

```bash
# 탐색·서술 분석 (반복 불필요)
cp configs/template/config_explore.yml configs/config_<실험명>.yml

# 또는 비교·통계 분석 (조건당 생물학적 샘플 ≥2)
cp configs/template/config_compare.yml configs/config_<실험명>.yml
```

| 목적 | 시작 템플릿 | 핵심 설정 |
|---|---|---|
| 탐색·서술 | [`configs/template/config_explore.yml`](configs/template/config_explore.yml) | `modules.*: false`, `optional_modules.exploration: true` |
| 비교·통계 | [`configs/template/config_compare.yml`](configs/template/config_compare.yml) | `modules.*: true`, `de_analysis.reference_level` |

> `configs/template/config.yml`은 Snakefile이 항상 자동으로 읽는 **기본값**입니다. 위 두 템플릿은 그 위에 덮어쓸 값만 담고 있으므로, 복사한 파일에서 필요한 항목만 수정하면 됩니다.

`configs/config_<실험명>.yml`을 열어 아래 항목을 수정합니다.

| 키 | 설명 |
|---|---|
| `single_cell.input_path` | 입력 파일 경로 (`.rds` 또는 `.h5ad`) |
| `single_cell.input_format` | `seurat` 또는 `anndata` |
| `single_cell.cell_type_column` | 세포 유형 컬럼명 |
| `single_cell.sample_column` | 샘플 ID 컬럼명 |
| `single_cell.condition_column` | 비교 조건 컬럼명 |
| `de_analysis.method` | `DESeq2` 또는 `edgeR` |
| `species` | `human` 또는 `mouse` |
| `modules.pseudobulk_de` | 비교 통계 사용 시 `true`, 탐색 전용이면 `false` |
| `modules.composition` | 세포 비율 비교 시 `true` |
| `modules.enrichment` | GO/KEGG 사용 시 `true` (`pseudobulk_de` 필요) |
| `optional_modules.exploration` | 유전자 패널 탐색 시 `true` |
| `optional_modules.trajectory` | `true` / `false` |
| `optional_modules.cell_communication` | `true` / `false` |

> `configs/config_*.yml`은 `.gitignore`로 제외되어 실험별 경로·조건 정보가 저장소에 올라가지 않습니다.

### 4. 파이프라인 실행

**비교 통계 포함 전체 실행 (Phase 1–6):**

```bash
snakemake --profile profiles/default \
  --configfile configs/config_<실험명>.yml \
  --cores 8
```

**탐색 전용 실행** (생물학적 반복 불필요, 탐색 템플릿을 복사한 경우):

```bash
snakemake --profile profiles/default \
  --configfile configs/config_<실험명>.yml \
  --cores 4
```

**Downstream만 재실행** (비교 분석에서 입력 로딩 생략, cutoff 변경 등):

```bash
snakemake --profile profiles/default \
  --configfile configs/config_<실험명>.yml \
  --config start_from=loaded \
  --cores 8
```

`profiles/default/config.yaml`에 `cores: 4`, `rerun-incomplete: true`,
`printshellcmds: true`가 정의되어 있습니다. `--cores`를 생략하면 4개 코어를 사용합니다.

> 기존 `run_pipeline.sh` / `run_downstream_only.sh`는 deprecated shim으로 남아 있지만,
> 새 실행에는 위 Snakemake 명령을 사용하세요.

---

## Config 주요 옵션

전체 기본값은 [`configs/template/config.yml`](configs/template/config.yml)에 있습니다. 모드별 시작점은 [`config_explore.yml`](configs/template/config_explore.yml), [`config_compare.yml`](configs/template/config_compare.yml)을 사용하세요.

```yaml
species: "human"          # "human" | "mouse"
output_dir: "output/my_experiment"   # 항상 output/<실험명>
start_from: "input"       # "input" | "loaded" (downstream-only)

padj_cutoff:   0.05       # DE 유의성 기준 (bulk 파이프라인과 동일 키)
log2fc_cutoff: 1.0

single_cell:
  input_path:      "path/to/seurat.rds"
  input_format:    "seurat"       # "seurat" | "anndata"
  cell_type_column: "cell_type"
  sample_column:   "sample_id"
  condition_column: "condition"
  cell_types:      null           # null = 전체 세포 유형
  min_cells_per_group: 10         # pseudobulk 최소 세포 수

de_analysis:
  method:          "DESeq2"       # "DESeq2" | "edgeR"
  reference_level: ""             # 빈 문자열 = 알파벳 첫 번째
  covariates:      []             # 예: ["batch", "sex"]

enrichment:
  go_ontologies:   ["BP", "CC", "MF"]
  run_kegg:        true
  pvalue_cutoff:   0.05
  qvalue_cutoff:   0.2

modules:
  pseudobulk_de: true    # 탐색 전용이면 false
  composition: true      # 탐색 전용이면 false
  enrichment: true       # 탐색 전용이면 false (pseudobulk_de 필요)

optional_modules:
  exploration:        false   # 탐색 전용이면 true
  trajectory:         false
  cell_communication: true

trajectory:
  method:         "monocle3"      # "monocle3" | "slingshot"
  root_cell_type: ""              # 필수 — 생물학적 판단 필요

cell_communication:
  method:               "CellChat"   # "CellChat" | "LIANA"
  conditions_to_compare: []          # 빈 배열 = 전체 조건

report:
  title:                 "scRNA-seq Analysis Report"
  author:                ""
  generate_methods_draft: true
```

---

## 출력 디렉토리 구조

`output_dir`은 **`output/<실험명>`** 형태로 지정합니다. 실행별 결과가 최상위 `output/` 아래에 서브폴더로 분리됩니다.

```
output/
└── <실험명>/                        # config의 output_dir (예: output/tutorial_compare)
    ├── config_used.yml              # 실행 시 config 스냅샷
    ├── loaded/
    │   └── seurat_validated.rds     # 검증된 Seurat 객체
    ├── pseudobulk_de/
    │   ├── de_<cell_type>.csv       # 세포 유형별 DE 결과
    │   ├── de_all_celltypes.csv     # 전체 통합 DE 결과
    │   ├── de_significant.csv       # 유의 유전자만 (padj + log2FC 필터)
    │   └── plots/
    │       └── volcano_<ct>_<contrast>.png
    ├── composition/
    │   ├── cell_type_proportions.csv
    │   ├── propeller_results.csv
    │   └── plots/
    │       ├── composition_bar.png
    │       └── composition_by_celltype.png
    ├── enrichment/
    │   ├── <ct>_<contrast>_GO_BP_up.csv
    │   ├── <ct>_<contrast>_KEGG_down.csv
    │   └── plots/
    │       └── <ct>/
    ├── trajectory/                  # optional_modules.trajectory: true 시
    │   ├── pseudotime_monocle3.csv
    │   ├── monocle3_cds.rds
    │   └── plots/
    ├── cell_communication/          # optional_modules.cell_communication: true 시
    │   ├── lr_interactions_<condition>.csv
    │   ├── cellchat_<condition>.rds
    │   └── plots/
    ├── tables/                      # canonical CSV table layer
    │   ├── manifest.json
    │   ├── index.tsv
    │   ├── README.md
    │   ├── analysis/datasets/*.csv
    │   └── expression/datasets/*.csv
    ├── export/
    │   ├── analysis/                # viewer analysis staging
    │   ├── expression/              # expression parquet + CSV source
    │   ├── excel/                   # dataset_type별 workbook
    │   │   ├── differential_expression.xlsx
    │   │   ├── go_analysis.xlsx
    │   │   ├── composition.xlsx
    │   │   ├── multi_group.xlsx
    │   │   ├── expression_summary.xlsx
    │   │   ├── excel_manifest.json
    │   │   └── README.md
    │   └── seqviewer/               # cmg-seqviewer Import Folder
    │       ├── seqviewer_manifest.json
    │       ├── datasets/*.parquet
    │       ├── export_index.tsv
    │       ├── COLUMNS.md
    │       └── export_complete.flag
    └── report/
        ├── report.html
        └── methods_draft.md         # 논문 Methods 초안 자동 생성
```

---

## 예시 데이터로 테스트

```bash
# Seurat 패키지에 내장된 pbmc3k 데이터 활용 예시
Rscript -e "
  library(Seurat)
  data('pbmc3k.final')
  obj <- pbmc3k.final
  obj[['sample_id']]  <- sample(c('S1','S2','S3','S4'), ncol(obj), replace=TRUE)
  obj[['condition']]  <- ifelse(obj[['sample_id', drop=TRUE]] %in% c('S1','S2'), 'ctrl', 'treat')
  obj[['cell_type']]  <- Idents(obj)
  saveRDS(obj, 'test_pbmc3k.rds')
"

cp configs/template/config_compare.yml configs/config_test.yml
# config_test.yml에서 input_path: "test_pbmc3k.rds" 로 수정 후:
snakemake --profile profiles/default --configfile configs/config_test.yml --cores 4
```

---

## Trajectory 분석 주의사항

Trajectory 분석은 **생물학적 판단이 필수**입니다.

```yaml
trajectory:
  root_cell_type: "Naive CD4 T"   # 반드시 명시 — 빈 문자열이면 오류
  root_cluster:   ""              # 같은 cell type 내 특정 cluster로 좁힐 때 사용
```

- `root_cell_type`은 분화의 시작점이 되는 세포 유형을 생물학적 지식에 근거하여 지정합니다.
- root를 잘못 지정하면 pseudotime 방향이 반전될 수 있습니다.
- root가 불명확한 경우에는 `optional_modules.trajectory: false`로 두고 별도 분석을 권장합니다.

---

## 방법론 참고

**Pseudobulk DE 권장 이유:**  
단일세포 수준에서 조건 간 DE를 직접 수행하면 세포 수가 많을수록 통계적 검정력이 과도하게 높아져 생물학적으로 의미 없는 유전자도 유의하게 나올 수 있습니다 (pseudo-replication). Pseudobulk 방식은 샘플을 통계적 단위로 사용하여 이 문제를 해결합니다.

**GeneRatio 해석:**  
clusterProfiler의 `GeneRatio`는 입력 유전자 중 해당 term에 속하는 비율을 나타내며, 유의성 기준이 아닙니다. 유의성은 `p.adjust` (padj)로 판단하세요.

---

## 소프트웨어 환경

| 소프트웨어 | 버전 | 용도 |
|---|---|---|
| R | ≥ 4.3 | 전체 분석 |
| Snakemake | ≥ 7.32 | 워크플로우 |
| Seurat | ≥ 5.0 | 단일세포 객체 |
| DESeq2 | Bioconductor | Pseudobulk DE |
| edgeR / limma | Bioconductor | Pseudobulk DE (alternative) |
| speckle (propeller) | Bioconductor | Compositional analysis |
| clusterProfiler | Bioconductor | GO / KEGG enrichment |
| CellChat | CRAN | Cell-cell communication |
| monocle3 | Bioconductor | Trajectory analysis |
| slingshot | Bioconductor | Trajectory analysis (alternative) |

---

## 주요 참고문헌

- Stuart T, et al. Comprehensive Integration of Single-Cell Data. *Cell* 2019.
- Love MI, et al. DESeq2. *Genome Biology* 2014.
- Robinson MD, et al. edgeR. *Bioinformatics* 2010.
- Phipson B, et al. propeller. *Nature Methods* 2022.
- Wu T, et al. clusterProfiler 4.0. *Innovation* 2021.
- Jin S, et al. CellChat. *Nature Communications* 2021.
- Cao J, et al. Monocle3. *Nature* 2019.
- Street K, et al. Slingshot. *BMC Genomics* 2018.
