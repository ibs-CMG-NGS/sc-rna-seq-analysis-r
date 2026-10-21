# HypoMap으로 시작하는 scRNA-seq 3차 분석

이 문서는 **생물학 전공자이지만 리눅스 터미널과 R 패키지 설치가 낯선 분**을 위한 실습 튜토리얼입니다.
공개 HypoMap atlas를 내려받아

1. 세포타입별 유전자 발현을 탐색하고,
2. 한 데이터셋 안에서 식이 조건(Normal chow vs Fasted)을 비교하는 통계 분석

까지 직접 실행합니다.

- 예상 소요: 환경 설치 30–60분, 데이터 다운로드 10–60분(네트워크 의존), 분석 실행 15분
- 권장 장비: Linux 또는 WSL2, **RAM 32GB 이상**, 디스크 여유 10GB 이상
- 검증 환경: R 4.5.3, Seurat 5.5.1, DESeq2 1.50.2, speckle 1.10.0, Snakemake 9.27.0
- 모든 명령은 저장소 최상위 디렉토리에서 실행합니다.

> 버전이 조금 다르면 출력 파일의 세부 수치가 달라질 수 있습니다. 명령과 흐름은 동일합니다.

---

## 0. 전체 흐름

```text
[1] conda 환경 준비
      ↓
[2] HypoMap 다운로드 (약 4GB)
      ↓
[3] 탐색적 분석: 유전자 패널 × 세포타입
      ↓
[4] 비교 분석: Dowsett10xnuc, Normal chow vs Fasted
      ↓
[5] HTML 리포트 + Methods 초안
```

각 단계 끝에는 **체크포인트**가 있습니다. 다음 단계로 넘어가기 전에 꼭 확인하세요.

### 명령을 붙여 넣는 방법

- 회색 코드 블록 전체를 복사해 터미널에 붙여 넣고 Enter를 누릅니다.
- `$` 기호는 프롬프트를 뜻하며, 실제로 입력하지 않습니다.
- 경로는 현재 저장소 기준 상대 경로입니다.
- Windows 사용자는 **WSL2 Ubuntu 터미널**에서 실행하세요. `C:\` 아래(`/mnt/c/...`)가 아니라 WSL 홈(`~/`)에 저장소를 두는 것을 권장합니다.

---

## 1. 준비물과 기본 개념

### 1.1 필요한 것

| 항목 | 권장 사양 | 확인 명령 |
|---|---|---|
| OS | Linux, WSL2, macOS | `uname -a` |
| RAM | 32GB 이상 | `free -h` |
| 디스크 | 10GB 이상 여유 | `df -h .` |
| 인터넷 | HypoMap 4.07GiB 다운로드 가능 | — |
| 소프트웨어 | Miniforge(conda), Git | `conda --version`, `git --version` |

16GB 노트북에서는 HypoMap 전체 객체를 R로 읽는 단계에서 메모리가 부족할 수 있습니다. 가능하면 분석 서버나 32GB 이상 장비를 사용하세요.

### 1.2 꼭 알아 둘 용어

| 용어 | 쉬운 설명 |
|---|---|
| 세포(cell) | 분석 단위. HypoMap에는 384,925개 세포가 있음 |
| 유전자(gene) | 발현을 측정하는 단위. HypoMap에는 57,362개 유전자가 있음 |
| count | 세포에서 특정 유전자 mRNA가 몇 개 읽혔는지 센 **원본 정수값** |
| UMI | 분자 하나를 세는 표식. count의 기본 단위 |
| Seurat object | 유전자×세포 행렬, 세포 정보(metadata), 차원축소 결과를 담은 R 객체 |
| metadata | 각 세포의 데이터셋, 샘플, 식이, 세포타입 같은 정보표 |
| cell type | 세포의 종류(예: Agrp 뉴런, 성상세포) |
| dataset | 원본 연구/실험 묶음. HypoMap은 18개 dataset을 통합 |
| pseudobulk | 같은 샘플·세포타입의 세포들을 count 합계로 묶은 것. bulk RNA-seq처럼 비교하기 위함 |
| log2FC | 두 조건의 발현 차이를 log2 배율로 표현. +1 = 2배 높음, -1 = 1/2 |
| padj / FDR | 수천 개 유전자를 동시에 검정할 때 우연을 보정한 p-value |
| composition | 조건에 따라 세포타입 **비율**이 달라지는지 보는 분석 |
| GO / KEGG | 유전자 집합이 어떤 생물학적 기능·경로에 몰려 있는지 보는 enrichment 분석 |

---

## 2. 환경 설치

### 2.1 터미널 열기

- **Linux**: `Ctrl + Alt + T`
- **macOS**: Spotlight에서 `Terminal` 검색
- **Windows**: PowerShell(관리자)에서 아래를 실행해 WSL2 Ubuntu를 설치한 뒤, 시작 메뉴에서 Ubuntu를 엽니다.

```powershell
wsl --install -d Ubuntu
```

설치 후 Ubuntu 터미널에서 나머지 명령을 실행합니다.

### 2.2 conda 라이선스 확인 (중요)

이 파이프라인은 **conda-forge + bioconda** 채널만 사용하고, Anaconda의 `defaults` 채널은 사용하지 않습니다.
`environment.yml`에도 `nodefaults`가 지정되어 있습니다.

먼저 현재 컴퓨터에 어떤 conda가 설치되어 있는지 확인합니다.

```bash
which conda
conda --version
conda info
conda config --show channels
```

`conda config --show channels` 결과에 `defaults`가 보이면 Anaconda 공식 저장소를 사용 중입니다.
Anaconda의 `defaults` 채널은 **조직 규모와 용도에 따라 유료 라이선스가 필요할 수 있습니다**(특히 200인 이상 조직).
소속 기관의 라이선스 정책을 먼저 확인하세요. 이 문서는 법률 자문이 아닙니다.

권장 선택:

| 현재 상태 | 권장 조치 |
|---|---|
| conda가 없음 | 아래 2.3에서 Miniforge 설치 |
| Miniforge 사용 중 (`conda info --base`가 `.../miniforge3`) | 그대로 사용 |
| Anaconda 사용 중, 기관 라이선스 보유 | 기존 conda 사용 가능. 단, `environment.yml`의 `nodefaults`는 유지 |
| Anaconda 사용 중, 라이선스 불명 | Miniforge를 별도 경로에 설치해 사용 |

`defaults`를 제거하려면 아래를 실행합니다. 이미 없으면 오류가 나므로 무시해도 됩니다.

```bash
conda config --remove channels defaults
```

채널 우선순위는 다음처럼 확인합니다.

```bash
conda config --show channel_priority
```

### 2.3 Miniforge 설치

Miniforge는 conda-forge 커뮤니티 채널을 기본으로 사용하는 conda 배포판입니다.

**Linux / WSL2**

```bash
curl -fsSL -o /tmp/Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash
exec bash -l
```

**macOS (Apple Silicon)**

```bash
curl -fsSL -o /tmp/Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh

bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init zsh
exec zsh -l
```

**macOS (Intel)**

`Miniforge3-MacOSX-arm64.sh` 대신 `Miniforge3-MacOSX-x86_64.sh`를 사용합니다.

설치 확인:

```bash
which conda
conda info --base
```

`conda info --base`가 `$HOME/miniforge3`를 가리키면 정상입니다.
기존 Anaconda가 먼저 잡히면 새 터미널을 열거나 `hash -r`을 실행하세요.

### 2.4 저장소 받기

```bash
git clone https://github.com/ibs-CMG-NGS/sc-rna-seq-analysis-r.git
cd sc-rna-seq-analysis-r
pwd
```

`pwd` 결과가 저장소 최상위 경로인지 확인하세요. 이후 모든 명령은 이 위치에서 실행합니다.

### 2.5 파이프라인 conda 환경 생성

```bash
conda env create -f environment.yml
conda activate scrnaseq-analysis
```

처음 설치에는 20–60분 정도 걸릴 수 있습니다. 중간에 멈춘 것처럼 보여도 기다려 주세요.

설치 확인:

```bash
Rscript -e 'for (p in c("Seurat","DESeq2","speckle","clusterProfiler","arrow")) cat(p, as.character(packageVersion(p)), "\n")'
snakemake --version
```

예상 출력(버전은 조금 다를 수 있음):

```text
Seurat 5.5.1
DESeq2 1.50.2
speckle 1.10.0
clusterProfiler 4.18.4
arrow 25.0.0
9.27.0
```

### 체크포인트 2

- [ ] `conda activate scrnaseq-analysis` 후 프롬프트 앞에 `(scrnaseq-analysis)`가 보인다
- [ ] `Rscript`와 `snakemake`가 실행된다
- [ ] Seurat, DESeq2, speckle이 로드된다

---

## 3. HypoMap 다운로드

### 3.1 HypoMap이란?

HypoMap은 생쥐 시상하부(hypothalamus)의 통합 단일세포 atlas입니다.

- 384,925 cells, 57,362 genes, 18 datasets
- scVI로 통합되어 있으며 `scvi`, `umap_scvi` 차원축소 제공
- 라이선스: CC BY 4.0
- DOI: 10.17863/CAM.87955
- 원 논문: Steuernagel et al., *Nature Metabolism* (2022)

### 3.2 다운로드

저장소 안에 `data/hypomap/`을 만들고 `curl`로 받습니다. `data/`는 `.gitignore`에 포함되어 있어 Git에 올라가지 않습니다.

```bash
mkdir -p data/hypomap
cd data/hypomap

curl -L --fail --retry 3 --retry-delay 5 -C - \
  -o hypoMap.rds \
  "https://www.repository.cam.ac.uk/bitstreams/904d34e8-7934-42fe-b637-edf0ac58c389/download"

cd ../..
```

- `-L`: 리다이렉트를 따라감
- `--fail`: HTTP 오류 시 실패 처리
- `-C -`: 중간에 끊겨도 이어받기
- 다운로드가 중단되면 **같은 명령을 다시 실행**하면 이어서 받습니다.

### 3.3 파일 크기와 무결성 확인

**Linux / WSL2**

```bash
ls -lh data/hypomap/hypoMap.rds
stat -c%s data/hypomap/hypoMap.rds
sha256sum data/hypomap/hypoMap.rds
```

**macOS**

```bash
ls -lh data/hypomap/hypoMap.rds
stat -f%z data/hypomap/hypoMap.rds
shasum -a 256 data/hypomap/hypoMap.rds
```

이 문서 작성 시점의 기준값:

```text
size   4369673011 bytes (4.07 GiB)
sha256 d0cee18fa726b5a20b16feb56689050700cd1dbc4fd6333df17a3753fcc8548d
```

크기가 다르면 다운로드가 덜 끝난 것입니다. 같은 `curl` 명령을 다시 실행하세요.
SHA256이 다르면 원본 저장소의 파일이 갱신되었을 수 있으니, 파일을 다시 받고 출처를 확인하세요.

### 3.4 R로 열리는지 확인

**RAM 22–26GB를 사용합니다.** 다른 무거운 프로그램을 닫고 실행하세요.

```bash
Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("data/hypomap/hypoMap.rds")
print(obj)
cat("cells:", ncol(obj), " genes:", nrow(obj), "\n")
'
```

예상 출력(핵심 부분):

```text
An object of class Seurat
57362 features across 384925 samples within 1 assay
Active assay: RNA (57362 features, 0 variable features)
 2 dimensional reductions calculated: scvi, umap_scvi
cells: 384925  genes: 57362
```

Seurat 5에서 v4 객체를 읽을 때 경고가 나올 수 있지만, 분석에는 문제가 없습니다.

### 체크포인트 3

- [ ] `data/hypomap/hypoMap.rds` 크기가 4,369,673,011 bytes
- [ ] SHA256이 위 값과 일치
- [ ] `readRDS` 결과가 57,362 × 384,925

---

## 4. 데이터 구조 이해

### 4.1 Seurat 객체 안을 들여다보기

```bash
Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("data/hypomap/hypoMap.rds")

cat("--- assays ---\n"); print(Assays(obj))
cat("--- reductions ---\n"); print(Reductions(obj))
cat("--- metadata columns ---\n"); print(colnames(obj@meta.data))
cat("--- datasets ---\n"); print(table(obj$Dataset))
cat("--- diets ---\n"); print(table(obj$Diet))
'
```

### 4.2 counts와 data의 차이

- `counts`: 원본 UMI count (정수). **pseudobulk 합산은 반드시 이 값을 사용**합니다.
- `data`: 정규화 후 log 변환된 값. 발현 비교·시각화에 사용합니다.

파이프라인은 `counts`가 정수인지 검사합니다. 정규화된 값을 counts 자리에 넣으면 오류가 납니다.

```bash
Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("data/hypomap/hypoMap.rds")
print(GetAssayData(obj, assay = "RNA", layer = "counts")[1:4, 1:4])
'
```

### 4.3 주요 metadata 컬럼

| 컬럼 | 의미 | 이 튜토리얼에서의 용도 |
|---|---|---|
| `Cell_ID` | 세포 ID | 식별자 |
| `Dataset` | 원본 연구 | 교란 확인, 탐색 분석의 분해 기준 |
| `Sample_ID` | 생물학적 샘플 | pseudobulk 단위 |
| `Diet` | 식이 조건 | 비교 조건 |
| `C66_named` | 66개 세포타입 라벨 | 비교 분석 |
| `C185_named` | 185개 세포타입 라벨 | 탐색 분석 |
| `Author_Exclude` | 원저자 제외 플래그 | `yes` 세포 제거 |
| `inferred_sex`, `Age` | 성별·나이 | 필요 시 공변량 |

`C66_named` 라벨은 `C66-46: Agrp.GABA-4`처럼 **클러스터 번호 + 대표 marker/세포타입** 형식입니다.
콜론과 공백이 들어가므로 YAML/config에 적을 때는 반드시 따옴표로 감싸세요.

### 4.4 Author_Exclude 처리

```bash
Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("data/hypomap/hypoMap.rds")
cat("excluded yes:", sum(obj$Author_Exclude == "yes", na.rm = TRUE), "\n")
keep <- is.na(obj$Author_Exclude) | obj$Author_Exclude != "yes"
cat("kept:", sum(keep), "\n")
'
```

예상:

```text
excluded yes: 28469
kept: 356456
```

탐색 스크립트와 export 스크립트는 이 플래그를 자동으로 처리합니다.
비교 분석은 아래 6.2의 subset 스크립트에서 미리 제거합니다.

### 4.5 C66과 C185

| 라벨 | 세포타입 수 | 특징 | 추천 용도 |
|---|---:|---|---|
| `C66_named` | 66 | 더 굵은 분류, 해석이 상대적으로 쉬움 | 조건 비교(DE) |
| `C185_named` | 185 | 더 세분화된 분류 | atlas 탐색 |

### 4.6 왜 전체 atlas로 DE를 돌리지 않나?

`Diet`는 `Dataset`과 강하게 얽혀 있습니다.

```bash
Rscript -e '
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("data/hypomap/hypoMap.rds")
print(table(obj$Dataset, obj$Diet))
'
```

대부분의 dataset은 한 가지 식이만 가지고 있습니다. 이런 상태에서 전체 atlas로 "Fasted vs Normal chow"를 비교하면
**식이 효과와 dataset/기술적 차이를 구분할 수 없습니다.**

그래서 이 튜토리얼은 **Dowsett10xnuc 한 dataset 안에서** 비교합니다.

| 조건 | 샘플 수 | 세포 수(제외 전) |
|---|---:|---:|
| Normal chow (adlib) | 4 | 24,955 |
| Fasted | 2 | 10,990 |

Fasted 샘플이 2개뿐이라 통계적 검정력이 낮습니다. 결과는 **가설 생성용**으로 해석해야 합니다.

### 체크포인트 4

- [ ] `counts`가 정수 행렬임을 확인
- [ ] `Author_Exclude` 처리 후 356,456 cells
- [ ] Dowsett10xnuc이 4 vs 2 샘플 구조임을 확인

---

## 5. 탐색적 분석: 유전자 패널 × 세포타입

### 5.1 목적

"이 유전자가 어느 세포타입에서 발현되는가?"를 **반복(replicate) 없이** 요약합니다.
이 단계는 p-value를 계산하지 않는 **서술적 분석**입니다.

예시 유전자 패널은 시상하부 에너지 대사와 관련된 유전자입니다.

```text
Glp1r, Pnoc, Agrp, Pomc, Mc4r, Hcrt, Oxt, Avp, Sst, Ghrh
```

### 5.2 설정 파일

`configs/tutorial/hypomap_explore.yml`을 사용합니다. 핵심 부분만 보면:

```yaml
output_dir: "output/tutorial_explore"

single_cell:
  input_path: "data/hypomap/hypoMap.rds"
  cell_type_column: "C185_named"
  sample_column: "Sample_ID"
  condition_column: "Diet"

modules:
  pseudobulk_de: false
  composition: false
  enrichment: false

optional_modules:
  exploration: true

exploration:
  enabled: true
  genes: ["Glp1r", "Pnoc", "Agrp", "Pomc", "Mc4r", "Hcrt", "Oxt", "Avp", "Sst", "Ghrh"]
  group_by: "C185_named"
  split_by: "Dataset"
  min_cells: 10
  exclude_column: "Author_Exclude"
  exclude_values: ["yes"]
```

- `modules`를 모두 `false`로 두어 반복이 필요한 분석은 건너뜁니다.
- `split_by: Dataset`은 같은 세포타입이라도 dataset별로 나누어 보여 줍니다.
- `min_cells: 10`은 10개 미만 세포의 그룹을 제외합니다.

### 5.3 실행

```bash
bash run_pipeline.sh configs/tutorial/hypomap_explore.yml --cores 4
```

예상 소요: **1–2분**, 최대 메모리 약 **21GB**.
(전체 HypoMap을 읽는 시간이 대부분입니다.)

### 5.4 출력

```text
output/
└── tutorial_explore/
    ├── config_used.yml
    ├── exploration/
    │   ├── gene_expression_by_celltype.csv
    │   ├── gene_expression_by_celltype_dataset.csv
    │   ├── top_celltypes.csv
    │   ├── summary.txt
    │   └── plots/
    │       ├── bar_Agrp.png
    │       ├── bar_Glp1r.png
    │       ├── ...
    │       ├── dotplot.png
    │       └── heatmap_top_celltypes.png
    └── report/
        └── report.html
```

### 5.5 결과표 읽기

`gene_expression_by_celltype.csv`의 컬럼:

| 컬럼 | 의미 |
|---|---|
| `gene` | 유전자 |
| `group` | 세포타입 |
| `n_cells` | 해당 세포타입의 세포 수 |
| `n_expressing` | 원본 count > 0인 세포 수 |
| `pct_expressing` | 발현 세포 비율(%) |
| `mean_expr` | 정규화 발현의 평균 |
| `median_expr` | 정규화 발현의 중앙값 |
| `max_expr` | 최대 발현 |
| `mean_counts` | 원본 count 평균 |
| `rank_by_mean` | 세포타입 간 평균 발현 순위 |
| `relative_to_max` | 최대 세포타입 대비 상대 발현 |

예시(Agrp):

| cell type | n_cells | pct_expressing | mean_expr |
|---|---:|---:|---:|
| `C185-115: Npy.Agrp.GABA-4` | 3,580 | 94.8% | 3.44 |
| `C185-116: Lepr.Agrp.GABA-4` | 554 | 89.7% | 2.19 |

Agrp가 Npy.Agrp 뉴런에서 높다는 것은 알려진 생물학과 일치합니다.
Glp1r은 `C185-113: Nkx2-4.GABA-3`에서 30.2%로 가장 높게 나옵니다.

### 5.6 그림 읽기

- `bar_<gene>.png`: 유전자별 상위 세포타입 bar plot
- `dotplot.png`: 점 크기 = 발현 세포 비율, 색 = 평균 발현
- `heatmap_top_celltypes.png`: 유전자 × 상위 세포타입 heatmap

### 5.7 dataset별 분해가 중요한 이유

`gene_expression_by_celltype_dataset.csv`는 `cell_type`과 `split`(Dataset)을 함께 보여 줍니다.
한 세포타입의 높은 발현이 **특정 dataset 한 곳에서만** 나오는지 확인할 수 있습니다.

예를 들어 Agrp는 `C185-115: Npy.Agrp.GABA-4`에서 여러 dataset에 걸쳐 높게 나타납니다.
반대로 어떤 유전자가 한 dataset에서만 높다면, 생물학적 차이보다 기술적 차이일 가능성을 먼저 의심해야 합니다.

### 5.8 해석할 때 주의

- 이 단계는 **p-value가 없습니다**. "차이가 있다/없다"를 통계적으로 말하지 않습니다.
- atlas는 여러 연구를 통합한 것이므로, dataset 간 발현 수준을 직접 비교하면 안 됩니다.
- 세포타입 라벨은 자동/수동 주석의 결과물이며 완벽하지 않습니다.
- `mean_expr`는 log-normalized 값이라 "2배" 같은 해석은 log2FC와 다릅니다.

### 체크포인트 5

- [ ] `output/tutorial_explore/exploration/`에 CSV 4개와 plot 12개 생성
- [ ] `summary.txt`에 `Cells used: 356456`, `185 groups` 표시
- [ ] `report/report.html`이 열림

---

## 6. 비교 통계 분석: Normal chow vs Fasted

### 6.1 설계

- 데이터: `Dowsett10xnuc` 단일 dataset
- 조건: `Normal chow` (reference) vs `Fasted`
- 샘플: 4 vs 2
- 세포타입: 해석이 쉽고 모든 샘플에서 10개 이상 세포가 있는 C66 타입 10개
- 분석: pseudobulk DESeq2 + propeller composition + GO/KEGG enrichment

> Fasted 샘플이 2개뿐이므로 결과는 **탐색적·가설 생성용**입니다. 논문의 확정 결론으로 사용하지 마세요.

### 6.2 비교용 subset 만들기

전체 atlas(4GB)를 매번 읽지 않도록, Dowsett 부분만 떼어 작은 파일로 저장합니다.

```bash
Rscript examples/hypomap/subset_dowsett.R
```

이 스크립트가 하는 일:

1. `Dataset == "Dowsett10xnuc"` 세포만 선택
2. `Author_Exclude == "yes"` 제거
3. 파이프라인이 기대하는 표준 컬럼 생성
   - `cell_type` ← `C66_named`
   - `sample_id` ← `Sample_ID`
   - `condition` ← `Diet` (levels: `Normal chow`, `Fasted`)
4. **모든 샘플에서 10개 이상** 세포가 있는 세포타입만 유지
5. `data/hypomap/dowsett_subset.rds`로 저장

예상 출력:

```text
[subset] Dowsett10xnuc: 35945 cells after exclusion
[subset] cell types kept: 43 / 63
Cells      : 33838
Genes      : 57362
Cell types : 43
Samples    : 6
Conditions :
Normal chow      Fasted
      24234        9604
```

- 소요: 약 2분
- 최대 메모리: 약 24GB (전체 atlas를 읽는 단계)
- 결과 파일: 약 494MB

### 6.3 설정 파일

`configs/tutorial/hypomap_dowsett_compare.yml`을 사용합니다. 핵심 부분:

```yaml
output_dir: "output/tutorial_compare"

single_cell:
  input_path: "data/hypomap/dowsett_subset.rds"
  cell_type_column: "cell_type"
  sample_column: "sample_id"
  condition_column: "condition"
  cell_types:
    - "C66-46: Agrp.GABA-4"
    - "C66-19: Pomc.GLU-5"
    - "C66-45: Ghrh.GABA-3"
    - "C66-4: Trh.GLU-2"
    - "C66-37: Crh.GABA-1"
    - "C66-47: Sst.GABA-4"
    - "C66-43: Nts.GABA-1"
    - "C66-28: Mixed.GABA-2"
    - "C66-53: Etnppl.Astrocytes"
    - "C66-55: Ermn.Oligodendrocytes"

modules:
  pseudobulk_de: true
  composition: true
  enrichment: true

de_analysis:
  method: "DESeq2"
  reference_level: "Normal chow"
```

세포타입 이름에 콜론과 공백이 있으므로 **각 항목을 따옴표로 감싸야** 합니다.

### 6.4 실행

```bash
bash run_pipeline.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4
```

예상 소요: **약 10분** (구성 요소별 실측)

| 단계 | 시간 |
|---|---:|
| 입력 로드·검증 | 24초 |
| composition (propeller) | 11초 |
| pseudobulk DE (DESeq2) | 43초 |
| GO/KEGG enrichment | 약 9분 |
| HTML 리포트 | 1초 |
| **전체** | **약 10분 20초** |

Enrichment가 대부분의 시간을 차지합니다. 정상입니다.

### 6.5 출력 구조

```text
output/
└── tutorial_compare/
    ├── config_used.yml
    ├── loaded/
    │   └── seurat_validated.rds
    ├── pseudobulk_de/
    │   ├── pseudobulk_counts.rds
    │   ├── pseudobulk_counts.tsv
    │   ├── de_all_celltypes.csv
    │   ├── de_significant.csv
    │   ├── de_<cell_type>.csv
    │   └── plots/
    │       └── volcano_<cell_type>_<contrast>.png
    ├── composition/
    │   ├── cell_type_proportions.csv
    │   ├── propeller_results.csv
    │   └── plots/
    │       ├── composition_bar.png
    │       └── composition_by_celltype.png
    ├── enrichment/
    │   ├── <cell_type>_GO_BP_up.csv
    │   ├── <cell_type>_GO_BP_down.csv
    │   ├── <cell_type>_KEGG_up.csv
    │   ├── <cell_type>_KEGG_down.csv
    │   └── plots/
    ├── export/
    │   └── seqviewer/               # cmg-seqviewer Import Folder
    │       ├── seqviewer_manifest.json
    │       ├── datasets/*.parquet
    │       ├── export_index.tsv
    │       ├── COLUMNS.md
    │       └── export_complete.flag
    └── report/
        ├── report.html
        └── methods_draft.md
```

### 6.6 Pseudobulk DE 결과 읽기

`pseudobulk_de/de_all_celltypes.csv`의 컬럼:

| 컬럼 | 의미 |
|---|---|
| `gene` | 유전자 |
| `baseMean` | 전체 평균 발현(count 기반) |
| `log2FoldChange` | Fasted / Normal chow의 log2 배율 |
| `lfcSE` | log2FC의 표준오차 |
| `stat` | 검정 통계량 |
| `pvalue` | 보정 전 p-value |
| `padj` | Benjamini-Hochberg 보정 p-value |
| `contrast` | 비교 이름 (`condition_Fasted_vs_Normal.chow`) |
| `cell_type` | 세포타입 |

해석:

- `log2FoldChange > 0`: **Fasted에서 더 높음**
- `log2FoldChange < 0`: **Normal chow에서 더 높음**
- `padj < 0.05` 그리고 `|log2FC| >= 1`을 만족하는 유전자를 significant로 저장

이 실행에서는:

- 전체 DE 행: 236,906
- significant 행: 1,441
- 세포타입별 significant 수:

| cell type | significant genes |
|---|---:|
| Agrp.GABA-4 | 384 |
| Sst.GABA-4 | 201 |
| Pomc.GLU-5 | 175 |
| Astrocytes | 129 |
| Oligodendrocytes | 120 |
| Ghrh.GABA-3 | 115 |
| Nts.GABA-1 | 97 |
| Trh.GLU-2 | 91 |
| Mixed.GABA-2 | 82 |
| Crh.GABA-1 | 47 |

예시로 Agrp 뉴런에서 padj 기준 상위 유전자는 `Sorcs1` (+2.51), `Pde10a` (+2.12), `Lepr` (+2.05)입니다.
Pomc 뉴런에서는 `Zbtb16` (+2.73), `Gal` (-8.05) 등이 보입니다.
**개별 유전자의 생물학적 의미는 별도 문헌 검토가 필요합니다.** 이 튜토리얼의 목적은 결과를 읽는 방법을 익히는 것입니다.

Volcano plot은 `pseudobulk_de/plots/`에 세포타입별로 저장됩니다.
x축은 log2FC, y축은 -log10(padj), 점선은 cutoff입니다.

### 6.7 Composition 결과 읽기

`composition/propeller_results.csv`:

| 컬럼 | 의미 |
|---|---|
| `cell_type` | 세포타입 |
| `PropMean.Normal.chow` | Normal chow에서의 평균 비율 |
| `PropMean.Fasted` | Fasted에서의 평균 비율 |
| `PropRatio` | 두 비율의 비 (Normal / Fasted) |
| `Tstatistic` | 검정 통계량 |
| `P.Value` | p-value |
| `FDR` | 보정 p-value |

이 실행에서는 선택한 10개 타입 중 **4개**가 FDR < 0.05였습니다.

| cell type | Normal chow | Fasted | FDR | 해석 |
|---|---:|---:|---:|---|
| Crh.GABA-1 | 6.39% | 2.26% | 0.025 | Normal chow에서 높음 |
| Nts.GABA-1 | 8.70% | 2.98% | 0.025 | Normal chow에서 높음 |
| Trh.GLU-2 | 17.86% | 5.59% | 0.035 | Normal chow에서 높음 |
| Ermn.Oligodendrocytes | 5.99% | 32.55% | 0.044 | Fasted에서 높음 |

주의: 이 비율은 **선택한 10개 세포타입 안에서의 상대 비율**입니다.
전체 43개 타입 기준 비율을 보려면 7.3의 "전체 타입 실행"을 참고하세요.

### 6.8 Enrichment 결과 읽기

파일 이름 규칙:

```text
<cell_type>_GO_<BP|CC|MF>_<up|down>.csv
<cell_type>_KEGG_<up|down>.csv
```

- `up`: Fasted에서 증가한 유전자
- `down`: Fasted에서 감소한 유전자
- `qvalue < 0.2`를 significant로 간주 (config에서 변경 가능)

모든 조합에 결과가 있는 것은 아닙니다. 유전자가 충분히 없거나 term이 없으면 빈 파일이 생성됩니다.
`enrichment/plots/`에는 세포타입별 dotplot이 저장됩니다.

Enrichment는 **DE 유전자 목록에 의존**합니다. DE cutoff를 바꾸면 enrichment 결과도 달라집니다.

### 6.9 HTML 리포트와 Methods 초안

`report/report.html`을 열면 다음 섹션이 순서대로 나옵니다.

1. Cell Type Composition
2. Pseudobulk Differential Expression
3. Enrichment Analysis
4. Session Info
5. Appendix: Methods Draft

`report/methods_draft.md`는 논문 Methods 초안입니다. 자동 생성된 초안이므로
실험 설계, 세포타입 주석 출처, 통계적 한계를 반드시 직접 검토·수정하세요.

### 6.10 해석할 때 주의

- **Fasted n=2**: 통계적 검정력이 낮고, 우연한 결과가 섞일 수 있습니다.
- **pseudobulk**: 세포 수가 많은 타입이 더 안정적입니다. 세포 수가 적은 타입은 결과가 흔들립니다.
- **composition**: 선택한 세포타입 집합에 따라 비율이 달라집니다.
- **enrichment**: 사용한 유전자 집합(GO/KEGG)과 background에 따라 결과가 달라집니다.
- **인과 해석 금지**: 이 분석은 조건 간 연관성을 보여 줄 뿐, 원인을 증명하지 않습니다.
- **dataset 특이성**: Dowsett10xnuc 한 연구의 결과입니다. 다른 dataset에서 재현되는지 확인해야 합니다.

### 체크포인트 6

- [ ] `de_all_celltypes.csv` 236,906행
- [ ] `de_significant.csv` 1,441행
- [ ] `propeller_results.csv` 10행, FDR < 0.05 4개
- [ ] `report/report.html`과 `methods_draft.md` 생성

---

## 7. 재실행과 파라미터 변경

### 7.1 입력 로딩을 건너뛰고 downstream만 다시 실행

cutoff만 바꿔 DE/enrichment를 다시 돌릴 때 사용합니다.

```bash
bash run_downstream_only.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4
```

이 스크립트는 `output/tutorial_compare/loaded/seurat_validated.rds`가 이미 있을 때만 동작합니다.

### 7.2 Cutoff 변경

config에서 수정합니다.

```yaml
padj_cutoff: 0.05
log2fc_cutoff: 1.0
```

더 엄격하게 하려면 `padj_cutoff: 0.01`, `log2fc_cutoff: 1.5`처럼 바꿉니다.
변경 후 `run_downstream_only.sh`로 재실행합니다.

### 7.3 전체 43개 세포타입으로 실행

본문은 해석을 쉽게 하기 위해 10개만 사용했습니다. 전체 robust 타입(43개)으로 돌리려면
config에서 `cell_types`를 `null`로 바꿉니다.

```yaml
single_cell:
  cell_types: null
```

그다음 실행합니다.

```bash
bash run_pipeline.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4
```

- DE 결과와 volcano plot이 43개로 늘어납니다.
- Enrichment 시간이 크게 증가합니다(수십 분).
- 리포트에는 volcano plot이 최대 20개까지만 삽입됩니다.

### 7.4 Reference 조건 바꾸기

```yaml
de_analysis:
  reference_level: "Fasted"
```

이 경우 `log2FoldChange > 0`은 **Normal chow에서 더 높음**을 뜻합니다.
해석 방향이 반대가 되므로 결과를 읽을 때 주의하세요.

### 7.5 처음부터 다시 실행

Snakemake는 완료된 파일을 재사용합니다. 강제로 다시 돌리려면:

```bash
bash run_pipeline.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4 -F
```

또는 출력 디렉토리를 지우고 실행합니다.

```bash
rm -rf output/tutorial_compare
bash run_pipeline.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4
```

---

## 8. 문제 해결

| 증상 | 원인 | 해결 |
|---|---|---|
| `conda: command not found` | 셸 초기화 안 됨 | `source ~/miniforge3/etc/profile.d/conda.sh` 후 `conda init bash`, 새 터미널 |
| `CondaHTTPError`, `403`, `defaults` 관련 오류 | Anaconda 채널/라이선스 | `conda config --show channels` 확인, `defaults` 제거, Miniforge 사용 |
| `conda env create`가 오래 걸림 | 패키지 다운로드/빌드 | 정상. 네트워크와 디스크 확인 |
| `Killed` 또는 R 프로세스 종료 | 메모리 부족 | 32GB 이상 장비 사용, 다른 프로그램 종료, WSL 메모리 증설 |
| WSL2에서 메모리 부족 | WSL 기본 메모리 제한 | Windows 사용자 홈에 `.wslconfig` 작성 후 `wsl --shutdown` |
| `curl` 다운로드 중단 | 네트워크 끊김 | 같은 `curl -C -` 명령 재실행 |
| SHA256 불일치 | 다운로드 손상/원본 갱신 | 파일 삭제 후 재다운로드, 출처 확인 |
| `Rscript: command not found` | conda 환경 미활성화 | `conda activate scrnaseq-analysis` |
| `there is no package called ...` | 환경 생성 불완전 | `conda env update -f environment.yml --prune` |
| `No cells found` | 세포타입 이름 오타/따옴표 누락 | config의 `cell_types`를 정확히 quoting, subset 파일의 라벨 확인 |
| `RDS` 로드가 너무 느림 | 4GB 객체 압축 해제 | 정상. SSD 사용 권장 |
| 리포트 그림이 안 보임 | 리포트만 따로 이동 | `report.html`은 self-contained이지만, 재생성 권장 |
| Enrichment가 오래 걸림 | GO/KEGG 조합이 많음 | 정상. `cell_types`를 줄이거나 `run_kegg: false` |
| `Permission denied` (`/mnt/c`) | Windows 파일시스템 권한 | WSL 홈(`~/`)에서 작업 |

### WSL2 메모리 설정 예시

Windows 사용자 홈(`C:\Users\<사용자>\.wslconfig`)에 작성:

```ini
[wsl2]
memory=28GB
processors=8
```

저장 후 PowerShell에서:

```powershell
wsl --shutdown
```

다시 Ubuntu를 엽니다. 호스트 RAM이 32GB라면 28GB 정도가 적당합니다.

---

## 9. 부록

### A. 용어집

| 용어 | 설명 |
|---|---|
| atlas | 여러 연구/데이터셋을 통합한 참조 지도 |
| batch effect | 실험·기술 차이로 생기는 인위적 패턴 |
| confounding | 관심 효과와 다른 요인이 뒤섞여 구분할 수 없는 상태 |
| DESeq2 | count 기반 차등발현 분석 패키지 |
| propeller | 세포 조성 비율을 검정하는 speckle 패키지의 함수 |
| pseudobulk | 세포들을 샘플 단위로 합쳐 만든 bulk 유사 발현 행렬 |
| ORA | 미리 정한 유전자 집합에 대한 over-representation analysis |
| FDR | false discovery rate, 다중검정 보정 오류율 |
| reference level | 비교의 기준이 되는 조건 |
| volcano plot | log2FC와 유의성을 함께 보여 주는 산점도 |

### B. 주요 출력 파일 요약

| 파일 | 내용 |
|---|---|
| `exploration/gene_expression_by_celltype.csv` | 유전자 × 세포타입 요약 |
| `exploration/gene_expression_by_celltype_dataset.csv` | 유전자 × 세포타입 × dataset 요약 |
| `pseudobulk_de/de_all_celltypes.csv` | 전체 DE 결과 |
| `pseudobulk_de/de_significant.csv` | cutoff 통과 DE 유전자 |
| `pseudobulk_de/pseudobulk_counts.rds` | pseudobulk count 행렬 |
| `composition/propeller_results.csv` | 세포 조성 검정 결과 |
| `enrichment/*.csv` | GO/KEGG enrichment 결과 |
| `report/report.html` | 통합 리포트 |
| `report/methods_draft.md` | Methods 초안 |

### C. cmg-seqviewer로 결과 보기 (선택)

이 파이프라인은 결과를 [cmg-seqviewer](https://github.com/ibs-CMG-NGS/cmg-seqviewer)에서 볼 수 있는
`seqviewer_manifest.json` + `datasets/*.parquet` 번들로 export할 수 있습니다.

예제 config(`configs/tutorial/hypomap_dowsett_compare.yml`)는 이미
`export_seqviewer.enabled: true`로 설정되어 있습니다. 다른 config에서 export하려면:

```yaml
export_seqviewer:
  enabled: true
```

그다음 다시 실행합니다.

```bash
bash run_pipeline.sh configs/tutorial/hypomap_dowsett_compare.yml --cores 4
```

출력:

```text
output/tutorial_compare/export/seqviewer/
├── seqviewer_manifest.json
├── datasets/
│   ├── DE__*.parquet
│   ├── ENR__*.parquet
│   ├── COMP__propeller.parquet
│   └── MG__*.parquet
├── export_index.tsv
├── COLUMNS.md
└── export_complete.flag
```

뷰어에서 이 폴더를 **Import Folder**로 선택합니다.

전체 유전자 × 세포타입 expression summary 번들을 만들려면 `expression_export` 모듈을 사용합니다.
`expression_export.enabled: true`로 설정하고 실행한 뒤:

```bash
Rscript scripts/08_build_viewer_bundle.R \
  --export_dir output/tutorial_compare/export/expression \
  --bundle_dir output/tutorial_compare/export/viewer/hypomap_dowsett \
  --prefix hypomap_dowsett \
  --alias "HypoMap Dowsett mean expression" \
  --cell_type_level C66 \
  --pct_scale 0-100
```

자세한 규격은 `README.md`의 "cmg-seqviewer 연동" 섹션을 참고하세요.

### D. GEO accession과 nf_fetchgeo의 경계

이 파이프라인은 **로컬 `.rds` 또는 `.h5ad` 파일만 입력**으로 받습니다.
GEO accession으로 직접 다운로드하는 기능은 포함되어 있지 않습니다.

GEO 데이터 수집은 별도 도구인 `nf_fetchgeo`(Nextflow)에서 수행하는 것을 권장합니다.

```text
nf_fetchgeo  →  로컬 .rds  →  이 파이프라인 (single_cell.input_path)
```

- `nf_fetchgeo`의 출력 경로를 `single_cell.input_path`에 지정합니다.
- 가능하면 GEO accession, 원본 파일명, SHA256, 다운로드 시각을 함께 보관하세요.
- `.rda`는 직접 지원하지 않습니다. 필요하면 R에서 읽어 `.rds`로 변환합니다.

```r
# .rda 파일을 받은 경우 (객체 이름은 파일마다 다를 수 있음)
load("GSEXXXXX_processed.rda")
print(ls())
saveRDS(obj, "GSEXXXXX_processed.rds")
```

`obj` 부분은 `ls()`로 확인한 실제 객체 이름으로 바꾸세요.
Snakefile 안에서 Nextflow를 호출하지 않고, **수집 단계와 분석 단계를 분리**하는 것이 재현성과 유지보수에 유리합니다.

### E. 재현성 체크리스트

- [ ] 사용한 conda 환경 이름과 `environment.yml` 버전 기록
- [ ] HypoMap 파일 크기와 SHA256 기록
- [ ] `config_used.yml` 보관 (파이프라인이 자동 저장)
- [ ] 분석 실행 날짜와 명령 기록
- [ ] `report/methods_draft.md`를 실제 분석에 맞게 수정
- [ ] 결과를 다른 dataset에서 재현할 수 있는지 확인
- [ ] `data/`, `output_*/`는 Git에 커밋하지 않기

### F. HypoMap 인용 및 라이선스

HypoMap을 사용한 결과를 발표할 때는 원 논문과 데이터 DOI를 인용하세요.

- Steuernagel L, et al. HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. *Nature Metabolism* (2022).
- DOI: 10.17863/CAM.87955
- 라이선스: CC BY 4.0

---

## 10. 다음 단계

- 다른 dataset에서 같은 비교를 해 보려면 `subset_dowsett.R`의 옵션을 바꿉니다.

  ```bash
  Rscript examples/hypomap/subset_dowsett.R \
    --dataset OtherDataset \
    --condition_column Treatment \
    --condition_levels "Control,Treated" \
    --output data/hypomap/other_subset.rds
  ```

  단, dataset × condition 교란과 샘플 수를 먼저 확인하고, 비교 config의
  `input_path`, `condition_column`, `reference_level`, `cell_types`도 새 데이터에 맞게 수정해야 합니다.
- 세포타입을 더 세분화하려면 `C185_named`를 사용합니다.
- 다른 유전자 패널을 보려면 `configs/tutorial/hypomap_explore.yml`의 `exploration.genes`를 수정합니다.
- 논문용 Methods는 `methods_draft.md`를 출발점으로 삼되, 반드시 직접 검토하세요.
