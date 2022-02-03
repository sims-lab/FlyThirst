#### 1. Load packages

    library(Seurat)
    library(SingleCellExperiment)
    library(ggplot2)
    library(Matrix)
    library(cowplot)
    library(tidyverse)
    library(BiocParallel)
    library(SummarizedExperiment)
    library(zinbwave)
    library(foreach)
    library(DESeq2)
    library(edgeR)
    library(eulerr)
    library(ggrastr)
    library(ggrepel)
    library(bookdown)
    setwd("~/Documents/Post_Doc/Drop-Seq/10X/Thirst2")

#### 2. Load data

    load("~/Documents/Post_Doc/Drop-Seq/10X/Thirst2/Thirst2_SCT_trimPlus.Robj")
    load("~/Documents/Post_Doc/Drop-Seq/10X/Thirst2/Thirst2_celltypes.Robj")

#### 3. Parameters

    pval <- 0.05
    logfc <- 1
    min.pct <- 0.25
    col_sat00 <- "#4e85c5"
    col_dep06 <- "#7e93a2"
    col_dep12 <- "#bb9e69"
    col_reh00 <- "#b373a6"

#### 4. ZINB-WaVE

For each cluster, calculate the observational weights according to the
ZINB model. There is a filtering step here, as only clusters with &gt;5
cells from each sample will be used. We later include these smaller
clusters into larger groups of cells to calculate differential
expression.  
This is time consuming so to be run on a HPC.

    mcp <- MulticoreParam(8)
    register(mcp)
    Thirst2_ZINB <- list()
    for(label in levels(Thirst2_SCT_trimPlus)){
      # Only run the analysis on clusters with at least 5 cells from each sample
      if(length(table(Thirst2_SCT_trimPlus[,Idents(Thirst2_SCT_trimPlus)==label]$orig.ident))==8 & min(table(Thirst2_SCT_trimPlus[,Idents(Thirst2_SCT_trimPlus)==label]$orig.ident))>=5){
        sce <- as.SingleCellExperiment(subset(Thirst2_SCT_trimPlus, cells=WhichCells(Thirst2_SCT_trimPlus, idents=label)))
        # only keep genes expressed in at least 10 cells 
        sce <- sce[rowSums(counts(sce) >= 1) >= 10,]
        # making sure the 'counts' assay is first, as we need to work with raw counts here
        nms <- c("counts", setdiff(assayNames(sce), "counts"))
        assays(sce) <- assays(sce)[nms]
        se <- SummarizedExperiment(assays=list(counts=as.matrix(counts(sce))), colData=colData(sce))
        # make design matrix
        se$stim <- factor(se$stim, levels=unique(se$stim))
        design <- model.matrix(~se$sex1*se$stim)
        # Calculate ZINB model.
        zinb <- zinbFit(Y=se, X=design, epsilon=1e6, K=2, verbose=T, BPPARAM=MulticoreParam(8))
        se_zinb <- zinbwave(Y=se, X=design, epsilon=1e6, fitted_model=zinb, K=2, verbose=T, observationalWeights=T)
        Thirst2_ZINB[[label]] <- se_zinb
      }
    }

#### 5. edgeR

Run edgeR on each cluster, using weights and latent factors calculated
with ZINB-WaVE.

    ZINB_edgeR <- list()
    for(i in names(Thirst2_ZINB)){
      exp <- Thirst2_ZINB[[i]]
      # Filter out genes expressed in less than 25% of the cells in each condition. 
      keep <- lapply(unique(exp$stim), function(x){counts(exp)[,grepl(x, colnames(counts(exp)))]})
      keep <- lapply(keep, function(x){rowSums(x!=0)/ncol(x)}) %>% do.call(cbind, .) %>% rowMax()>min.pct
      exp <- exp[keep,]
      # Make dataset & design
      dge <- DGEList(assay(exp))
      dge <- calcNormFactors(dge)
      design <- model.matrix(~sex1+stim, data=colData(exp))
      dge$weights <- assay(exp, "weights")
      # Estimate dispersion. 
      dge <- estimateDisp(dge, design)
      # Because this is single-cell data, it is better to use the normal likelihood ratio test glmFit.
      ZINB_edgeR[[i]] <- glmFit(dge, design)
    }

For each gene, calculate logFC and p-value between 12h-dehydrated and
sated flies (dep12 vs sat00). Genes with abs(log2FC) &gt; 1 and adjusted
p-value &lt; 0.05 are considered differently expressed.

    ZINB_edgeR.table.ALL <- data.frame()
    for(i in names(ZINB_edgeR)){
      lrt <- glmWeightedF(ZINB_edgeR[[i]], coef = 4)
      res <- topTags(lrt, 1e6) 
      res.t <- res$table %>% rownames_to_column("Gene") %>% mutate(Cluster=i, comparison=res$comparison, threshold=ifelse(abs(logFC)>1 & FDR<.05, "padj < 0.05 &\nlog2FC > 1", ifelse(abs(logFC)>1, "log2FC > 1", ifelse(FDR<.05, "padj < 0.05", "Not different")))) %>% mutate(nsig=sum(grepl("&", threshold)), nhits=sum(.$threshold!="Not different"))
      ZINB_edgeR.table.ALL <- rbind(ZINB_edgeR.table.ALL, res.t)
    }
    save(ZINB_edgeR.table.ALL, file="~/Documents/Post_Doc/Drop-Seq/10X/Thirst2/ZINB_edgeR.table.ALL.Robj")

#### 6. DESeq2

Run DESeq2 on each cluster, using weights and latent factors calculated
with ZINB-WaVE.

    ZINB_DESeq <- list()
    for(i in names(Thirst2_ZINB)){
      cat(i, "\n")
      exp <- Thirst2_ZINB[[i]]
      # Filter out genes expressed in less than 25% of cells in each conditions. 
      keep <- lapply(unique(exp$stim), function(x){counts(exp)[,grepl(x, colnames(counts(exp)))]})
      keep <- lapply(keep, function(x){rowSums(x!=0)/ncol(x)}) %>% do.call(cbind, .) %>% rowMax()>min.pct
      exp <- exp[keep,]
      # Run DESeq
      dds <- DESeqDataSet(exp, design=~sex1+stim)
      scr <- scran::computeSumFactors(dds)
      sizeFactors(dds) <- sizeFactors(scr)
      dds <- DESeq(dds, test="LRT", reduced=~sex1, minmu=1e-6, minReplicatesForReplace=Inf)
      ZINB_DESeq[[i]] <- dds
    }

For each gene, calculate logFC and p-value between 12h-dehydrated and
sated flies (dep12 vs sat00) using lfcshrink. Genes with abs(log2FC)
&gt; 1 and adjusted p-value &lt; 0.05 are considered differently
expressed.

    ZINB_DESeq.table.ALL <- data.frame()
    for(i in names(ZINB_DESeq)){
      res <- lfcShrink(ZINB_DESeq[[i]], contrast=c("stim", "dep12", "sat00"), type = "normal", alpha=.05)
      res$padj <- p.adjust(res$pvalue, method="BH")
      res.t <- res %>% as.data.frame() %>% rownames_to_column("Gene") %>% mutate(Cluster=i, comparison=gsub(" ", "_", sub("^.*stim ", "", res@elementMetadata$description[2])), threshold=ifelse(abs(log2FoldChange)>1 & padj<.05, "padj < 0.05 &\nlogFC > 1", ifelse(abs(log2FoldChange)>1, "logFC > 1", ifelse(padj<.05, "padj < 0.05", "Not different")))) %>% mutate(nsig=sum(grepl("&", threshold)), nhits=sum(.$threshold!="Not different"))
      ZINB_DESeq.table.ALL <- rbind(ZINB_DESeq.table.ALL, res.t)
    }
    save(ZINB_DESeq.table.ALL, file="~/Documents/Post_Doc/Drop-Seq/10X/Thirst2/ZINB_DESeq.table.ALL.Robj")

#### 7. Run the analysis on groups of clusters

Measure differential expression in groups of relevant clusters.

    Thirst2_regrouped <- list()
    IDs <- levels(Thirst2_SCT_trimPlus)
    Thirst2_regrouped[["abKCs"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("abKCs", IDs)])]
    Thirst2_regrouped[["prKCs"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("primeKCs", IDs)])]
    # Gamma KCs are already all grouped together
    Thirst2_regrouped[["PAM"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("PAM", IDs)])]
    Thirst2_regrouped[["DANs_other"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("DA$", IDs)])]
    Thirst2_regrouped[["TA"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("TA", IDs)])]
    # OA neurons are already all grouped together.
    # Serotonergic neurons are already all grouped together
    Thirst2_regrouped[["CortexGlia"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("Cortex", IDs)])]
    Thirst2_regrouped[["EnsheathingGlia"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("Ensheath", IDs)])]
    Thirst2_regrouped[["SurfaceGlia"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("Surface", IDs)])]
    Thirst2_regrouped[["Astro"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("Astro", IDs)])]
    Thirst2_regrouped[["FBEB"]] <- Thirst2_SCT_trimPlus[,WhichCells(Thirst2_SCT_trimPlus, idents=IDs[grepl("(FB|EB)", IDs)])]

##### ZINB-WaVE

    cl <- makeCluster(30)
    registerDoParallel(cl)
    Thirst2_regrouped_ZINB <- foreach(label=names(Thirst2_regrouped)) %dopar% {
      if(length(table(Thirst2_regrouped[[label]]$orig.ident))==8 & min(table(Thirst2_regrouped[[label]]$orig.ident))>=5){
        sce <- as.SingleCellExperiment(Thirst2_regrouped[[label]])
        sce <- sce[rowSums(counts(sce) >= 1) >= 10,]
        nms <- c("counts", setdiff(assayNames(sce), "counts"))
        assays(sce) <- assays(sce)[nms]
        se <- SummarizedExperiment(assays=list(counts=as.matrix(counts(sce))), colData=colData(sce))
        se$stim <- factor(se$stim, levels=unique(se$stim))
        design <- model.matrix(~se$sex1*se$stim)
        zinb <- zinbFit(Y=se, X=design, epsilon=1e6, K=2, verbose=T, BPPARAM=MulticoreParam(3))
        se_zinb <- zinbwave(Y=se, X=design, epsilon=1e6, fitted_model=zinb, K=2, verbose=T, observationalWeights=T)
      }
    }
    names(Thirst2_regrouped_ZINB) <- names(Thirst2_regrouped)

##### edgeR

    ZINB_regr_edgeR <- list()
    for(i in names(Thirst2_regrouped_ZINB)){
      cat(i, "\n")
      exp <- Thirst2_regrouped_ZINB[[i]]
      keep <- lapply(unique(exp$stim), function(x){counts(exp)[,grepl(x, colnames(counts(exp)))]})
      keep <- lapply(keep, function(x){rowSums(x!=0)/ncol(x)}) %>% do.call(cbind, .) %>% rowMax()>min.pct
      exp <- exp[keep,]
      dge <- DGEList(assay(exp))
      dge <- calcNormFactors(dge)
      design <- model.matrix(~sex1+stim, data=colData(exp))
      dge$weights <- assay(exp, "weights")
      dge <- estimateDisp(dge, design)
      ZINB_regr_edgeR[[i]] <- glmFit(dge, design)
    }
    ZINB_edgeR.table.regrouped <- data.frame()
    for(i in names(ZINB_regr_edgeR)){
      lrt <- glmWeightedF(ZINB_regr_edgeR[[i]], coef=4)
      res <- topTags(lrt, 1e6) 
      res.t <- res$table %>% rownames_to_column("Gene") %>% mutate(Cluster=i, comparison=res$comparison, threshold=ifelse(abs(logFC)>1 & FDR<.05, "padj < 0.05 &\nlog2FC > 1", ifelse(abs(logFC)>1, "log2FC > 1", ifelse(FDR<.05, "padj < 0.05", "Not different")))) %>% mutate(nsig=sum(grepl("&", threshold)), nhits=sum(.$threshold!="Not different"))
      ZINB_edgeR.table.regrouped <- rbind(ZINB_edgeR.table.regrouped, res.t)
    }

##### DESeq2

    ZINB_DESeq.table.regrouped <- data.frame()
    for(i in names(Thirst2_regrouped_ZINB)){
      exp <- Thirst2_regrouped_ZINB[[i]]
      keep <- lapply(unique(exp$stim), function(x){counts(exp)[,grepl(x, colnames(counts(exp)))]})
      keep <- lapply(keep, function(x){rowSums(x!=0)/ncol(x)}) %>% do.call(cbind, .) %>% rowMax()>min.pct
      exp <- exp[keep,]
      dds <- DESeqDataSet(exp, design=~sex1+stim)
      scr <- scran::computeSumFactors(dds)
      sizeFactors(dds) <- sizeFactors(scr)
      dds <- DESeq(dds, test="LRT", reduced=~sex1, minmu=1e-6, minReplicatesForReplace=Inf)
      res <- lfcShrink(dds, contrast=c("stim", "dep12", "sat00"), type = "normal")
      res$padj <- p.adjust(res$pvalue, method="BH")
      res.t <- res %>% as.data.frame() %>% rownames_to_column("Gene") %>% mutate(Cluster=i, comparison=gsub(" ", "_", sub("^.*stim ", "", res@elementMetadata$description[2])), threshold=ifelse(abs(log2FoldChange)>1 & padj<.05, "padj < 0.05 &\nlog2FC > 1", ifelse(abs(log2FoldChange)>1, "log2FC > 1", ifelse(padj<.05, "padj < 0.05", "Not different")))) %>% mutate(nsig=sum(grepl("&", threshold)), nhits=sum(.$threshold!="Not different"))
      ZINB_DESeq.table.regrouped <- rbind(ZINB_DESeq.table.regrouped, res.t)
    }
