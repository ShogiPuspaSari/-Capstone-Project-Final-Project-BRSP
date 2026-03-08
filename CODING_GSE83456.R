#Modul: Analisis Transcriptomics TB
#Dataset: GSE83456 (Tuberculosis)
#Platform: Microarray (Illumina GPL10558)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)

# PART A. PENGANTAR KONSEP
# Analisis ekspresi gen membandingkan tingkat ekspresi gen
# antara dua kondisi biologis (TB aktif vs Healthy Control)
# Metode yang digunakan adalah limma (Linear Models for Microarray Data)

# PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE) 
#Apa itu package? 
#Package adalah kumpulan fungsi siap pakai di R
#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor 

#1. Install BiocManager (manajer paket Bioconductor) 
#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”
# Install BiocManager

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Install paket Bioconductor (GEOquery & limma) 
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 


# Install annotation package sesuai platform
# GPL10558 = Illumina HumanHT-12 V4.0
BiocManager::install("illuminaHumanv4.db", ask = FALSE, update = FALSE)

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}


#4. Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
# Load library
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 


#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE83456", GSEMatrix = TRUE, AnnotGPL = FALSE)

if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))


#PART D. PRE-PROCESSING DATA EKSPRESI 
ex <- exprs(gset)

qx <- as.numeric(quantile(ex,
                          c(0.,0.25,0.5,0.75,0.99,1.0),
                          na.rm=TRUE))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}


# PART E. DEFINISI KELOMPOK SAMPEL
# contoh pembagian kelompok:
# TB aktif vs Healthy Control

gsms <- "111111111111000000000000"
sml <- strsplit(gsms, split="")[[1]]

gs <- factor(sml)

groups <- make.names(c("TB","Healthy"))
levels(gs) <- groups
gset$group <- gs

# PART F. DESIGN MATRIX
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ]

v <- vooma(gset, design, plot = TRUE)
v$genes <- fData(gset)

# PART G. ANALISIS DIFFERENTIAL EXPRESSION
fit <- lmFit(v)

cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2,
               adjust="fdr",
               sort.by="B",
               number=250)

write.table(tT,
            file=stdout(),
            row.names=FALSE,
            sep="\t")


# Histogram P-value
tT2 <- topTable(fit2,
                adjust="fdr",
                sort.by="B",
                number=Inf)

hist(tT2$adj.P.Val,
     col="grey",
     border="white",
     xlab="P-adj",
     ylab="Number of genes",
     main="Adjusted P-value Distribution")

# DEG classification
dT <- decideTests(fit2,
                  adjust.method="fdr",
                  p.value=0.05,
                  lfc=0)

vennDiagram(dT, circle.col=palette())


# QQ plot
t.good <- which(!is.na(fit2$F))

qqt(fit2$t[t.good],
    fit2$df.total[t.good],
    main="Moderated t statistic")

# Volcano plot
ct <- 1

volcanoplot(fit2,
            coef=ct,
            main=colnames(fit2)[ct],
            pch=20,
            highlight=length(which(dT[,ct]!=0)),
            names=rep('+', nrow(fit2)))

# MD plot
plotMD(fit2,
       column=ct,
       status=dT[,ct],
       legend=FALSE,
       pch=20,
       cex=1)

abline(h=0)

# PART H. VISUALISASI DATA
ex <- exprs(gset)

ord <- order(gs)

palette(c("#1B9E77","#7570B3"))

par(mar=c(7,4,2,1))

title <- paste("GSE83456","/",annotation(gset))

boxplot(ex[,ord],
        boxwex=0.6,
        notch=TRUE,
        main=title,
        outline=FALSE,
        las=2,
        col=gs[ord])

legend("topleft",
       groups,
       fill=palette(),
       bty="n")


# Density plot
par(mar=c(4,4,2,1))

title <- paste("GSE83456","/",annotation(gset),"value distribution")

plotDensities(ex,
              group=gs,
              main=title,
              legend="topright")

# UMAP plot
ex <- na.omit(ex)
ex <- ex[!duplicated(ex),]

ump <- umap(t(ex),
            n_neighbors=11,
            random_state=123)

par(mar=c(3,3,2,6), xpd=TRUE)

plot(ump$layout,
     main="UMAP plot",
     col=gs,
     pch=20,
     cex=1.5)

legend("topright",
       inset=c(-0.15,0),
       legend=levels(gs),
       pch=20,
       col=1:nlevels(gs),
       title="Group")

