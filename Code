
#### 1 remove inter-batch differences ####
if (!requireNamespace("sva", quietly = TRUE)) {  
  install.packages("sva")  
  library(sva)  
} else {  
  library(sva)  
}  
 
matrix1 <- read.delim("matrix1.txt", row.names=1, check.names=FALSE)  
matrix2 <- read.delim("matrix2.txt", row.names=1, check.names=FALSE)  
  
batch1 <- matrix1$Batch  
batch2 <- matrix2$Batch  

all_genes <- rbind(matrix1[, -which(names(matrix1) == "Batch"), drop = FALSE],  
                   matrix2[, -which(names(matrix2) == "Batch"), drop = FALSE])  
all_batches <- c(batch1, batch2 + max(batch1)) # 确保批次编号不重叠  

all_matrix <- data.frame(all_genes, Batch = as.factor(all_batches))  

adjusted_data <- combat(dat=as.matrix(all_genes), batch=all_batches)  

write.table(adjusted_data, "adjusted_matrix.txt", sep="\t", quote=FALSE, row.names=TRUE)  

adjusted_with_batch <- data.frame(adjusted_data, Batch=all_batches)  
write.table(adjusted_with_batch, "adjusted_matrix_with_batch.txt", sep="\t", quote=FALSE, row.names=TRUE)


###### 2 Identification of Differentially Expressed ERGs ####
if (!requireNamespace("limma", quietly = TRUE)) {  
  install.packages("limma")  
}  
library(limma)  
  
expr_matrix <- read.csv("expression_matrix.csv", row.names=1)  

group_info <- read.csv("group_info.csv")  
  
group_factor <- factor(group_info$Group, labels=c("con", "treat"))  

design <- model.matrix(~0 + group_factor)  
colnames(design) <- levels(group_factor)  

contrast.matrix <- makeContrasts(treat-con, levels=design)  

fit <- lmFit(expr_matrix, design)  

fit2 <- contrasts.fit(fit, contrast.matrix)  
fit2 <- eBayes(fit2)  

topTable <- topTable(fit2, coef=1, adjust.method="fdr", number=Inf, sort.by="B", p.value=0.05, lfc=0.5)  

write.csv(topTable, "DEGs.csv", row.names=FALSE)  

DE_genes <- rownames(topTable[topTable$adj.P.Val < 0.05 & abs(topTable$logFC) > 0.5, ])  
print(DE_genes)

##### 3 ssGSEA ######
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="merge without gsm114232.txt"         #表达输入文件
gmtFile="immune.gmt"            #免疫数据集文件
clusterFile="cluster.txt"    #分型结果文件
setwd("")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,"cluster",drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", fill="cluster",
            xlab="", 
            ylab="Immune infiltration",
            legend.title="ERGscluster",
            palette=bioCol)
p=p+rotate_x_text(50)
#输出图形文件
pdf(file="boxplot.pdf", width=6, height=6)
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()

##### 4 Consensus clustering ######
library(ConsensusClusterPlus)      #引用包
expFile="ERGsExp.txt"          #表达数据文件
workDir=""     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

#删掉对照组样品，只保留实验组样品
group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="treat"]

#聚类
maxK=5
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=10,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="pearson",
                             seed=123456,
                             plot="png")


#输出分型结果
clusterNum=2        #分几类，根据前面的图形判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)

##### 5 PCA #####
#引用包
library(limma)
library(ggplot2)

clusterFile="cluster.txt"     #分型的结果文件
setwd("")      #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
cluster=as.vector(rt[,ncol(rt)])

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], cluster=cluster)
PCA.mean=aggregate(PCA[,1:2], list(cluster=PCA$cluster), mean)

#设置颜色
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
CluCol=bioCol[1:length(levels(factor(cluster)))]


#定义椭圆函数
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$cluster))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$cluster==g,],
                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),
                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,
                                                                   center=c(mean(PC1),mean(PC2))))), cluster=g))
}

#绘制PCA图形
pdf(file="PCA.pdf", height=3, width=4)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = cluster)) +
  scale_colour_manual(name="cluster", values =CluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=cluster), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$cluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##### 6 Calculation of the Cluster Score #####
expFile="ERGsExp.txt"      #表达数据文件
setwd("")     #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

#PCA分析
pca=prcomp(data, scale=TRUE)
value=predict(pca)
score=value[,1]+value[,2]
score=as.data.frame(score)
scoreOut=rbind(id=colnames(score), score)
write.table(scoreOut, file="score.txt", sep="\t", quote=F, col.names=F)

##### 7 Nomogram #####
#引用包
library(rms)
library(rmda)

inputFile="input.txt"       #输入文件
setwd("")      #设置工作目录

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type~ Age+Blood_Sugar+Complications+Gender+NEUT+SIAH2, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of DF")
#输出列线图
pdf("Nom.pdf", width=6, height=4)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5, height=5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

#绘制决策曲线
rt$Type=ifelse(rt$Type=="con", 0, 1)
dc=decision_curve(Type ~ Age+Blood_Sugar+Complications+Gender+NEUT+SIAH2, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
#输出DCA图形
pdf(file="DCA.pdf", width=5, height=5)
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

#绘制临床影响曲线
pdf(file="clinical_impact.pdf", width=5, height=5)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("red", "blue"),
                     legend="bottomleft")
dev.off()


###### 8 Seurat #####

library(Seurat)  

counts <- read.csv("matrix.csv", row.names = 1) # 读取数据  

seurat_obj <- CreateSeuratObject(counts = counts)  

seurat_obj <- FilterCells(seurat_obj, min.cells = 3, min.features = 200) 
seurat_obj <- FilterCells(seurat_obj, idents = Idents(seurat_obj), percent.mito = 10) # 移除线粒体基因比例超过10%的细胞  

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")  

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", min.mean.scale = 0.0125)  

seurat_obj <- ScaleData(seurat_obj)  
seurat_obj <- RunPCA(seurat_obj, dims = 1:1500) # 使用前1500个高变基因进行PCA  

seurat_obj <- RunUMAP(seurat_obj, dims = 1:80, resolution = 0.5) # 使用前80个主成分进行UMAP降维  
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:80) # 使用前80个主成分进行邻域计算  
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # 基于邻域图进行聚类  



##### 9 monocle ####
mmt<-sqy_cell[,Idents(sqy_cell)%in%c("HS","LS")]

mmt$cell_type_val<-Idents(mmt)
mmt@meta.data[1:5,]

mmt<-FindVariableFeatures(mmt)

matrix<-as.matrix(mmt@assays$RNA@counts)  #提取最原始的count数据
dim(matrix) #看多少行多少列
matrix[1:500,1:10]

gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))

mmt@meta.data[1:5,]
sample_ann <- mmt@meta.data

fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)  #建立monocle研究对象，并进行归一化和质控
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)

sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 10)) #必须至少在10个细胞里（num_cells_expressed）有表达，才纳入该基因,可以自己调
fData(sc_cds_2)[1:5,]

ordering_genes<-mmt@assays$RNA@var.features

sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
#plot_ordering_genes(sc_cds2) 
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")   #降到2维
sc_cds2 <- orderCells(sc_cds2)  #把cell的顺序排出来
saveRDS(sc_cds2,"monocle_for_plot.rds")


##### 10 cellchat #####
sqy <- readRDS("")
sqy<-createCellChat(sqy)
sqy@DB <- CellChatDB.human
sqy<-subsetData(sqy)

sqy <- identifyOverExpressedGenes(sqy)  #确定每个细胞亚群中的过表达基因
sqy <- identifyOverExpressedInteractions(sqy) #寻找过表达的interaction
sqy <- projectData(sqy, PPI.human)  #向ppi投射
sqy<- computeCommunProb(sqy, raw.use = T,) #算interaction的可能性

sqy <- filterCommunication(sqy, min.cells = 10)  #去除interaction很少的细胞
sqy <- computeCommunProbPathway(sqy)  #计算通路
sqy <- netAnalysis_computeCentrality(sqy, slot.name = "netP") #信号网络的一个拓扑学分析

saveRDS(sqy,file="sqy_cellchat.rds")

#### This code is by Qingyuan Sun in Shandong University
#### Date : June 14th, 2024
