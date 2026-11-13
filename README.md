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



