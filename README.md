# scRNA-seq Downstream Analysis Pipeline

R과 Snakemake 기반의 재현 가능한 단일세포 RNA-seq 3차 분석 파이프라인입니다.  
Seurat (`.rds`) 또는 AnnData (`.h5ad`) 형식을 입력으로 받아, 세포 조성 분석 → Pseudobulk DE → 유전자 세트 풍부화 → (선택) Trajectory → (선택) 세포 간 커뮤니케이션 → 통합 HTML 보고서까지 자동 생성합니다.

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

## cmg-seqviewer 연동 (export)

`export_seqviewer` rule이 `output/export/seqviewer/`에 뷰어 번들을 만듭니다 (config의 `export_seqviewer.enabled`로 끌 수 있음).

- `seqviewer_manifest.json`: 데이터셋 목록 (`datasets` 배열, `dataset_type` 명시)
- `datasets/*.parquet`: 데이터셋 파일
- `export_index.tsv`: 파일별 유형·세포타입·조건 목록
- `COLUMNS.md`: 컬럼 규약 (schema_version 기록)
- `export_complete.flag`: 완료 표시

| 파일 | dataset_type | 내용 |
|---|---|---|
| `DE__<cell_type>__<contrast>.parquet` | differential_expression | 세포타입 × contrast별 DE 결과 |
| `ENR__<cell_type>_<contrast>__GO_<ONT>_<up\|down>.parquet`, KEGG 동일 패턴 | go_analysis | GO/KEGG 결과 (세포타입 × contrast × 방향) |
| `COMP__propeller.parquet` | composition | 세포타입 비율 검정 (propeller) |
| `MG__<cell_type>.parquet` | multi_group | pseudobulk raw count, 샘플 열은 `<condition>_<rep>` |

뷰어에서는 `output/export/seqviewer/` 폴더 자체를 **Import Folder**로 선택합니다. `metadata.json`이 아니라 `seqviewer_manifest.json`을 읽습니다.

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

출력: `output/label_transfer/` 아래 `annotated.rds` (predicted + prediction_probability 반영, `cell_type` 컬럼이 없었던 경우에만 predictions로 채움) + `label_provenance.csv` (세포별 라벨·확률) + `label_provenance.yml` (모델 경로, 패키지 버전, input 형태 등 재현성 정보).
</details>

**bulk RNA-seq 파이프라인과의 일관성 유지:**
- `configs/template/config.yml` (Git 추적) + `configs/config_*.yml` (gitignore)
- `padj_cutoff`, `log2fc_cutoff` 최상위 레벨 키
- `run_pipeline.sh` / `run_downstream_only.sh` 실행 스크립트
- `config_used.yml` 자동 저장, 출력물 CSV + PNG

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

```bash
cp configs/template/config.yml configs/config_<실험명>.yml
```

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
| `optional_modules.trajectory` | `true` / `false` |
| `optional_modules.cell_communication` | `true` / `false` |

> `configs/config_*.yml`은 `.gitignore`로 제외되어 실험별 경로·조건 정보가 저장소에 올라가지 않습니다.

### 4. 파이프라인 실행

**전체 실행 (Phase 1–6):**

```bash
bash run_pipeline.sh configs/config_<실험명>.yml --cores 8
```

**Downstream만 재실행** (입력 로딩 생략, cutoff 변경 등):

```bash
bash run_downstream_only.sh configs/config_<실험명>.yml --cores 8
```

**Snakemake 직접 실행:**

```bash
snakemake --configfile configs/config_<실험명>.yml --cores 8 --rerun-incomplete
```

---

## Config 주요 옵션

```yaml
species: "human"          # "human" | "mouse"
output_dir: "output"

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

optional_modules:
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

```
output/
├── config_used.yml                  # 실행 시 config 스냅샷
├── loaded/
│   └── seurat_validated.rds         # 검증된 Seurat 객체
├── pseudobulk_de/
│   ├── de_<cell_type>.csv           # 세포 유형별 DE 결과
│   ├── de_all_celltypes.csv         # 전체 통합 DE 결과
│   ├── de_significant.csv           # 유의 유전자만 (padj + log2FC 필터)
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
├── trajectory/                      # optional_modules.trajectory: true 시
│   ├── pseudotime_monocle3.csv
│   ├── monocle3_cds.rds
│   └── plots/
├── cell_communication/              # optional_modules.cell_communication: true 시
│   ├── lr_interactions_<condition>.csv
│   ├── cellchat_<condition>.rds
│   └── plots/
├── export/
│   └── seqviewer/                   # cmg-seqviewer Import Folder
│       ├── seqviewer_manifest.json
│       ├── datasets/*.parquet
│       ├── export_index.tsv
│       ├── COLUMNS.md
│       └── export_complete.flag
└── report/
    ├── report.html
    └── methods_draft.md             # 논문 Methods 초안 자동 생성
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

cp configs/template/config.yml configs/config_test.yml
# config_test.yml에서 input_path: "test_pbmc3k.rds" 로 수정 후:
bash run_pipeline.sh configs/config_test.yml --cores 4
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
