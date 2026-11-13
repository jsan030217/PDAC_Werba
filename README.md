# PDAC_Werba
pdac 11 patients werba sample with R


##1. 관련 패키지 loading
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ade4)
library(scCustomize)
library(RColorBrewer)
library(viridis)


##2. 데이터 로딩 및 메타데이터 추가 후 통합
#데이터: pdac_werba data(11 samples)
#2-1. 데이터 불러오기: 데이터 이름<-Read10X("데이터 주소 경로")
#지정 디렉토리 안에는 matrix.mtx(count값), barcode.tsv(세포 바코드), feature.tsv(유전자 이름)의 세가지 파일이 있으며, 이를 Read10X로 읽어 유전자 X 세포 형태의 희소행렬로 만들어 준다.
p4.data <- Read10X("C:/Users/user/Desktop/jsdata/P4")
p5.data <- Read10X("C:/Users/user/Desktop/jsdata/P5")
p7.data <- Read10X("C:/Users/user/Desktop/jsdata/P7")
p9.data <- Read10X("C:/Users/user/Desktop/jsdata/P9")
p13.data <- Read10X("C:/Users/user/Desktop/jsdata/P13")
p15.data <- Read10X("C:/Users/user/Desktop/jsdata/P15")
p19.data <- Read10X("C:/Users/user/Desktop/jsdata/P19")
p20.data <- Read10X("C:/Users/user/Desktop/jsdata/P20")
p22.data <- Read10X("C:/Users/user/Desktop/jsdata/P22")
p23.data <- Read10X("C:/Users/user/Desktop/jsdata/P23")
p26.data <- Read10X("C:/Users/user/Desktop/jsdata/P26")

#2-2. 1차 필터링 및 seurat으로 객체 변환: 데이터 이름<-createseuratobject(counts=읽을 데이터, project=객체의 프로젝트 이름, min.cells=3, min.features=100) 
#counts에는 어떤 데이터를 불러올지, project에는 어떤 이름으로 orig.ident metada에 넣을지, 최소 세포 및 최소 검출 유전자 수를 설정하여 이보다 아래의 세포 및 유전자는 1차적으로 filtering
p4 <- CreateSeuratObject(counts = p4.data, project = 'P4', min.cells = 3, min.features = 100)
p5 <- CreateSeuratObject(counts = p5.data, project = 'P5', min.cells = 3, min.features = 100)
p7 <- CreateSeuratObject(counts = p7.data, project = 'P7', min.cells = 3, min.features = 100)
p9 <- CreateSeuratObject(counts = p9.data, project = 'P9', min.cells = 3, min.features = 100)
p13 <- CreateSeuratObject(counts = p13.data, project = 'P13', min.cells = 3, min.features = 100)
p15 <- CreateSeuratObject(counts = p15.data, project = 'P15', min.cells = 3, min.features = 100)
p19 <- CreateSeuratObject(counts = p19.data, project = 'P19', min.cells = 3, min.features = 100)
p20 <- CreateSeuratObject(counts = p20.data, project = 'P20', min.cells = 3, min.features = 100)
p22 <- CreateSeuratObject(counts = p22.data, project = 'P22', min.cells = 3, min.features = 100)
p23 <- CreateSeuratObject(counts = p23.data, project = 'P23', min.cells = 3, min.features = 100)
p26 <- CreateSeuratObject(counts = p26.data, project = 'P26', min.cells = 3, min.features = 100)

#2-3. metadata 추가
#위에서 만들어진 P(n) data에 ID, patient_ID, disease_state 등 추가적인 샘플 정보 기입
p4$ID <- "PDAC41"
p5$ID <- "PDAC42"
p7$ID <- "PDAC43"
p9$ID <- "PDAC44"
p13$ID <- "PDAC45"
p15$ID <- "PDAC46"
p19$ID <- "PDAC47"
p20$ID <- "PDAC48"
p22$ID <- "PDAC49"
p23$ID <- "PDAC50"
p26$ID <- "PDAC51"

p4$patient_ID <- "werba_P4"
p5$patient_ID <- "werba_P5"
p7$patient_ID <- "werba_P7"
p9$patient_ID <- "werba_P9"
p13$patient_ID <- "werba_P13"
p15$patient_ID <- "werba_P15"
p19$patient_ID <- "werba_P19"
p20$patient_ID <- "werba_P20"
p22$patient_ID <- "werba_P22"
p23$patient_ID <- "werba_P23"
p26$patient_ID <- "werba_P26"

p4$disease_state <- "PDAC"
p5$disease_state <- "PDAC"
p7$disease_state <- "PDAC"
p9$disease_state <- "PDAC"
p13$disease_state <- "PDAC"
p15$disease_state <- "PDAC"
p19$disease_state <- "PDAC"
p20$disease_state <- "PDAC"
p22$disease_state <- "PDAC"
p23$disease_state <- "PDAC"
p26$disease_state <- "PDAC"

#2-4. merge 함수 통합
#데이터 이름<-merge(첫번째 객체, y = c(나머지 객체)) 
#merge는 첫번째 객체에 나머지 객체를 단순히 이어 붙여서 세포가 전부 모인 하나의 객체로 만들어지며, 유전자는 합집합으로 처리된다. 이는 단순 이어 붙이는 것이라 이후 환자 간 batch effect는 따로 처리해주어야 한다.(integration 함수는 통합 및 batch effect 한번에 해결)
pdac_werba<-merge(p4, y = c(p5, p7, p9, p13, p15, p19, p20, p22, p23, p26)) 

#1번째 세이브포인트
saveRDS(pdac_werba, file = "C:/Users/user/Desktop/jsdata/merged_raw.rds")


##3. doublet remove
#3-1. 관련 패키지 loading
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)

#3-2. merge로 인한 샘플별 레이어 분리 문제 발생->JoinLayers() 로해 
pdac_werba<-JoinLayers(pdac_werba)

#3-3. seurat을 sce 형태로 변환
pdac_werba.sce <- as.SingleCellExperiment(pdac_werba)

#3-4. scDblFinder 알고리즘으로 샘플별 patient_ID로 하여 doublet 확인
#이때 samples="patient_ID" 로 하여 같음 샘플 안에서의 더블렛으로 나눠서 분석함(batch effect 해결 전이기 때문에 따로 분석해야 함)
pdac_werba.sce <- scDblFinder(pdac_werba.sce, samples = "patient_ID")

#3-5. singlet 과 doublet 확인
table(pdac_werba.sce$scDblFinder.class)
<img width="959" height="1024" alt="1  remov doublet table" src="https://github.com/user-attachments/assets/22b81ddf-1d07-4945-b332-5e4ba8919a56" />

#3-6. sce를 다시 seurat 형태로 변환
#기존에는 data="logcounts" 였지만, 에러 발생(data에 logcounts값이 없어서)->data=NULL로 바꿔서 실행
pdac_werba_seurat <- as.Seurat(pdac_werba.sce, counts = "counts",
                               data = NULL,
                               verbose = TRUE)

#3-7. originalexp를 RNA로 되돌리는 과정-이는 Seurat::Assays(pdac_werba_seurat) 코드로 확인 결과 "RNA"로 나와서 생략함.
pdac_werba_seurat <- RenameAssays(object = pdac_werba_seurat, originalexp = "RNA")

#3-8. 실제 doublet 제거 과정(singlet으로 분류되는 것만 보존)
pdac_werba_2 <- subset(pdac_werba_seurat, subset = scDblFinder.class == "singlet")

#2번째 세이브포인트
saveRDS(pdac_werba_2, file = "C:/Users/user/Desktop/jsdata/merged_singlet.rds")


##4. HGNC symbol matching
#4-1. 관련 패키지 loading
library(HGNChelper)

#4-2. sample의 gene에 대해서 같지만 구버전이거나 다른 이름으로 된 gene을 찾는다.
#HGNChelper::checkGeneSymbols() 함수로 유전자 symbol 검증함수를 호출하고 입력된 유전자 목록을 공식 HGNC 데이터베이스와 대조한다.
#이때 불러오는 pdac_werba_2@assays$RNA@counts@Dimnames[[1]] 는 gene name이다.(왜 rownames(pdac_werba_2) 를 사용하지 않으시는지)
#옵션 중 unmapped.as.na=FALSE 는 매핑되지 않은 유전자를 NA(결측치)로 두지 않고 기존 유전자로 둔다는 의미이다.
HGNC.pdac_werba_2 <- HGNChelper::checkGeneSymbols(pdac_werba_2@assays$RNA@counts@Dimnames[[1]], unmapped.as.na = FALSE, map = NULL)

#4-3. HGNC symbol gene으로 교체
#RNA assays 부분만 임시 변수로 복사해 오고, counts와 data 행렬의 행 이름을 최신 표준 심볼로 교체하고, 이 RNA부분을 기존의 pdac_werba_2 객체에 다시 덮어씌워 저장한다.
RNA <- pdac_werba_2@assays$RNA
RNA@counts@Dimnames[[1]]            <- HGNC.pdac_werba_2$Suggested.Symbol
RNA@data@Dimnames[[1]]                <- HGNC.pdac_werba_2$Suggested.Symbol
pdac_werba_2@assays$RNA <- RNA

#4-4. 결과 확인
#만약 HGNC symbol로 모든 유전자가 정확히 매칭되었다면 2번째 코드에서 2개 이상 등장하는 유전자는 없으므로 FALSE가 나와야 한다.
table(pdac_werba_2@assays$RNA@counts@Dimnames[[1]])
any(table(pdac_werba_2@assays$RNA@counts@Dimnames[[1]]) > 1)
#결과가 TRUE로 나와서 추가 코드 시행(make.unique())
<img width="507" height="768" alt="symbol_matching 실패 gene 모음" src="https://github.com/user-attachments/assets/30df77fd-b44c-4f05-936b-fb4d630d5bfb" />
clean_symbols<-make.unique(HGNC.pdac_werba_2$Suggested.Symbol)
RNA<-pdac_werba_2@assays$RNA
RNA@counts@Dimnames[]<-clean_symbols
RNA@data@Dimnames[]<-clean_symbols
rownames(RNA@meta.features)<-clean_symbols
pdac_werba_2@assays$RNA<-RNA
#이후 FALSE로 나와서 다음 단계 진행


##5. Filtering
#nFeature, nCount, percent.mt를 보면서 QC
#5-1. percent.mt 계산
pdac_werba_2[["percent.mt"]] <- PercentageFeatureSet(pdac_werba_2, pattern = "^MT-")

#5-2. 시각화
VlnPlot(pdac_werba_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
[vlnplot.pdf](https://github.com/user-attachments/files/34789401/vlnplot.pdf)
plot1 <- FeatureScatter(pdac_werba_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pdac_werba_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
[P1+P2.pdf](https://github.com/user-attachments/files/34789411/P1%2BP2.pdf)

#5-3. plot을 토대로 filtering 실시
pdac_werba_2 <- subset(pdac_werba_2, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)

#5-4. filtering 후 plot
[vlnplot(filtering 후).pdf](https://github.com/user-attachments/files/34789442/vlnplot.filtering.pdf)
[P3+P4(filtering 후).pdf](https://github.com/user-attachments/files/34789448/P3%2BP4.filtering.pdf)


##6. processing(normalization, HVG selection, scaling, runPCA)
#normalization: sequencing depth 차이 보정하기 위해서. 세포마다 전체 RNA 함량이 다를 수 있기 때문에 
pdac_werba_2 <- NormalizeData(pdac_werba_2, normalization.method = "LogNormalize", scale.factor = 10000)
#Highly variable genes(HVG) selection: 설정한 nFeatures 수만큼 세포 간 차이가 뚜렷한 고변이 유전자를 선별한다.
pdac_werba_2 <- FindVariableFeatures(pdac_werba_2, selection.method = "vst", nfeatures = 3000)
#all.genes<-rownames(pdac_werba_2) 를 사용하여 추후 DoHeatmap() 함수로 마커 유전자들의 발현 패턴을 그릴 때, HVG 목록에 없던 유전자도 스케일링된 값으로 깔끔하게 시각화 가능
all.genes <- rownames(pdac_werba_2)
#scaling: 모든 유전자에 대해서 평균을 0으로, 분산을 1로 맞춘다. 이를 통해 발현량이 원래 높은 유전자가 PCA나 heatmap에서 결과를 독점하는 현상을 방지하고, 모든 유전자가 동등한 가중치를 갖도록 만든다.
pdac_werba_2 <- ScaleData(pdac_werba_2, features = all.genes)
#PCA: 고변이유전자 3000개에 대해서만 PCA를 수행한다. PCA의 결과 3000개의 유전자로 그린 3000차원의 고차원에서, 데이터의 분포에 따라 분산이 가장 큰 방향의 PC1을 잡고, 그와 직교하는 PC2 를 만들며 새로운 축을 만든다. 이렇게 PCA를 수행하여 나온 PC를 몇 개 사용할지 결정하여 umap을 그리며, 그때는 elbowplot을 이용하여 PC의 parameter을 결정할 수 있다. 이렇게 고차원에서 저차원으로 차원 축소를 하여 실제 분석에 용이하게 사용한다. 
pdac_werba_2 <- RunPCA(pdac_werba_2, features = VariableFeatures(object = pdac_werba_2))
<img width="484" height="671" alt="3  RunPCA 결과_PC의 positive negative_1" src="https://github.com/user-attachments/assets/0cbccad5-a772-4b6d-b1df-9b7ddaf4ae4b" />
<img width="487" height="424" alt="3  RunPCA 결과_PC의 positive negative_2" src="https://github.com/user-attachments/assets/de185c98-e82a-4891-8f5b-a0a07e5c672f" />


##7. Harmony & clustering 
#Harmony: 여러 환자 데이터들을 merge로 단순 병합하였기 때문에 아직 환자마다 개인차로 인한 batch effect를 제거한다.
#UMAP: 보정된 harmony공간의 상위 30개 차원을 사용해 세포들을 2차원 평면에 시각화할 좌표를 계산한다. (reduction='harmony' 로 지정하지 않으면 기본값인 pca가 나온다.)
#FindNeighbors: 보정된 30개 차원 공간에서 각 세포와 유전자 발현 패턴이 가장 비슷한 이웃 세포들을 찾아 KNN그래프로 연결
#FindClusters: 연결망을 바탕으로 세포들을 clustering(resolution의 수치가 낮을수록 큰 덩어리 위주로 굵직하게 나뉜다. 0.1 의 파라미터로 설정해 major cell types를 확인할 때 사용한다.)
library(harmony)
pdac_werba_2 <- RunHarmony(pdac_werba_2, group.by.vars = "patient_ID")
pdac_werba_2 <- RunUMAP(pdac_werba_2, reduction = "harmony", dims = 1:30)
pdac_werba_2 <- FindNeighbors(pdac_werba_2, reduction = "harmony", dims = 1:30)
pdac_werba_2 <- FindClusters(pdac_werba_2, resolution = 0.1)

#3번째 세이브 포인트
saveRDS(pdac_werba_2, file = "C:/Users/user/Desktop/jsdata/merged_werba_processed.rds")

#데이터 불러오기
setwd("C:/Users/user/Desktop/jsdata/result&plot/Werba_PDAC")
werba <- readRDS("merged_werba_processed.rds")
werba <- FindClusters(werba, resolution = 0.1)

#8. visualization














