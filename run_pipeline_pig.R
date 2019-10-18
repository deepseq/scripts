library(rgl)
library(dendextend)
library(scater)
library(limma)
library(scran)
library(matrixStats)
library(gplots)
library(scatterplot3d)
library(pvclust)
library(ggplot2)
library(easyGgplot2)
library(scde)
library(DESeq2)
library(goseq)

source("heatmap.2.nolayout.R", chdir = TRUE)
'%nin%' <- Negate('%in%')

#################################################################################################################################################
#                                                                   function                                                                    #
#################################################################################################################################################

random.matrix  <- matrix(runif(500, min = -1, max = 1), nrow = 50)
quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.01)
color.palette  <- colorRampPalette(c("#91CF60", "#FFFFBF", "#FC8D59"))(length(palette.breaks) - 1)

legend.col <- function(col, lev)
{
    opar <- par
    n <- length(col)
    bx <- par("usr")
    
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
    bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
    
    xx <- rep(box.cx, each = 2)
    
    par(xpd = TRUE)
    for(i in 1:n){
        
        yy <- c(box.cy[1] + (box.sy * (i - 1)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
        
    }
    par(new = TRUE)
    plot(0, 0, type = "n",
    ylim = c(min(lev), max(lev)),
    yaxt = "n", ylab = "",
    xaxt = "n", xlab = "",
    frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = .25)
    par <- opar
}

run.scde <- function(input, treatnum, controlnum, treatstr, controlstr, output_dir, output_prefix, maplot)
{
    dir.create(output_dir)
    
    cond <- c(rep(treatstr,treatnum), rep(controlstr,controlnum))
    names(cond) <- colnames(input)
    groups <- factor(cond, levels = c(treatstr,controlstr))
    
    ifm <- scde.error.models(counts = input, groups = groups, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
    
    valid.cells <- ifm$corr.a > 0
    ifm <- ifm[valid.cells, ]
    prior <- scde.expression.prior(models = ifm, counts = input, length.out = 400, show.plot = FALSE)
    
    ediff <- scde.expression.difference(ifm, input, prior, groups = groups, n.randomizations = 100, n.cores = 1, verbose = 1)
    
    p.values <- 2*pnorm(abs(ediff$Z), lower.tail = F) # 2-tailed p-value
    p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail = F) # Adjusted to control for FDR
    sig <- which(p.values.adj < 0.05)
    
    ord <- order(p.values.adj[sig]) # order by p-value
    de <- cbind(ediff[sig,1:3], p.values.adj[sig])[ord,]
    colnames(de) <- c("Lower bound","log2 fold change","Upper bound","p-value")
    
    out <- paste(output_dir,"/",output_prefix,".scde_sig.txt", sep = "")
    
    write.table(de, file = out, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    if(maplot == "T")
    {
        macol <- input[,1]
        names(macol) <- rownames(input)
        macol.check <- names(macol)%in%rownames(de)
        macol[macol.check==TRUE] <- "red"
        macol[macol.check==FALSE] <- "grey"
        plot(rowMeans(log2(input+1)), ediff[rownames(input),2], pch = 16, cex = 0.5, col = macol, xlab = "mean normailised expression", ylab = "log2 of fold change", ylim = c(-10,10))
        abline(h = c(-1,1), col = "dodgerblue", lwd = 2)
    }
    
    return(de)
}

run.goseq <- function(deg, genelen, goanno, geneidmap, output_dir, output_prefix)
{
    deg.check <- rownames(genelen)%in%rownames(deg)
    names(deg.check) <- rownames(genelen)
    pwf <- nullp(deg.check, bias.data = genelen[,1])
    wall <- goseq(pwf, gene2cat = goanno)
    enriched <- wall$category[wall$over_represented_pvalue < 0.05]
    deg.govis <- subset(wall, over_represented_pvalue < 0.05, select = c(category,over_represented_pvalue))
    deg.gores <- subset(wall, over_represented_pvalue < 0.05)
    
    deg.gogenes <- goanno[goanno[,1]%in%rownames(deg),]
    deg.gogenes <- deg.gogenes[order(deg.gogenes$V2),]
    
    deg.gogenes.enrich <- deg.gogenes[deg.gogenes[,2]%in%deg.gores[,1],]
    
    nrows <- dim(deg.gores)[1]
    ncols <- dim(deg.gores)[2]
    
    for(idx in 1:nrows)
    {
        tmp <- deg.gogenes.enrich[deg.gogenes.enrich[,2]%in%deg.gores[idx,1],1]
        tmp <- as.matrix(tmp)
        tmp <- geneidmap[tmp,2]
        tmp <- as.matrix(tmp)
        deg.gores[idx,ncols+1] <- paste(tmp, collapse = ",")
    }
    
    colnames(deg.gores)[ncols+1] <- "GOGenes"
    
    deg.govis.out <- paste(output_dir,"/",output_prefix,".scde_govis.txt", sep = "")
    deg.gores.out <- paste(output_dir,"/",output_prefix,".scde_gores.txt", sep = "")
    deg.gogenes.out <- paste(output_dir,"/",output_prefix,".scde_gogenes.txt", sep = "")
    
    write.table(deg.govis, file = deg.govis.out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(deg.gores, file = deg.gores.out, quote = FALSE, row.names = FALSE, sep = "\t")
    write.table(deg.gogenes.enrich, file = deg.gogenes.out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
    return(deg.gores)
}


#################################################################################################################################################
#                                                                   files                                                                       #
#################################################################################################################################################

pigcell_counts <- read.table("allSample.counts.geneName_ID.txt", header = T, row.names = 1)
pigcell_counts.value <- pigcell_counts[,c(3:dim(pigcell_counts)[2])]
pigcell.geneid <- pigcell_counts[,c(1,2)]

pigcell_counts_unigene <- read.table("allSample.unigene.samtoolscounts.geneName_ID.txt", header = T, row.names = 1)
pigcell_counts_unigene.value <- pigcell_counts_unigene[,c(3:dim(pigcell_counts_unigene)[2])]
pigcell_unigene.geneid <- pigcell_counts_unigene[,c(1,2)]

pigcell_counts.value <- rbind(pigcell_counts.value, pigcell_counts_unigene.value[,colnames(pigcell_counts.value)])
pigcell.geneid <- rbind(pigcell.geneid, pigcell_unigene.geneid)

pigcell_tpms <- read.table("allSample.tpms.geneName_ID.txt", header = T, row.names = 1)
pigcell_tpms.value <- pigcell_tpms[,c(3:dim(pigcell_tpms)[2])]

pigcell_tpms_unigene <- read.table("allSample.unigene.tpms.geneName_ID.txt", header = T, row.names = 1)
pigcell_tpms_unigene.value <- pigcell_tpms_unigene[,c(3:dim(pigcell_tpms_unigene)[2])]

pigcell_tpms.value <- rbind(pigcell_tpms.value, pigcell_tpms_unigene.value[,colnames(pigcell_tpms.value)])

pigcell_mtcounts <- read.table("allSample.mtcounts.geneName.txt", header = T, row.names = 1)
pigcell_mtcounts.value <- pigcell_mtcounts[,c(2:dim(pigcell_mtcounts)[2])]

pigcell_stagecol <- read.table("samplecolor.txt", row.names = 1)
pigcell_stagecol <- as.matrix(pigcell_stagecol)

pigcell_id <- read.table("sampleid_cellid.txt", header = T, row.names = 1)
pigcell_id <- as.matrix(pigcell_id)

pigcell_entrezid <- read.table("Sus_scrofa_ensembl87_v10_2_clean.geneIDvsEntrezID.format.txt", row.names = 1)
pigcell_entrezid <- as.matrix(pigcell_entrezid)

pigcell_entrezid_re <- read.table("Sus_scrofa_ensembl87_v10_2_clean.geneIDvsEntrezID.reformat.txt", row.names = 1)
pigcell_entrezid_re <- as.matrix(pigcell_entrezid_re)

pig_goterm <- read.table("Sus_scrofa_ensembl87_v10_2_clean.geneNameIDvsGO.txt", sep = '\t', comment = "")
pig_genelen <- read.table("Sus_scrofa_ensembl87_v10_2_clean.geneLength.geneNameID.txt", row.names = 1, sep = '\t', comment = "")


#################################################################################################################################################
#                                                           kick out unwanted cells                                                             #
#################################################################################################################################################

pigcell_unwanted <- c("i702_i503", "i710_i503", "i701_i504", "i707_i504")

pigcell_counts.value <- pigcell_counts.value[,!(colnames(pigcell_counts.value)%in%pigcell_unwanted)]
pigcell_tpms.value <- pigcell_tpms.value[,colnames(pigcell_counts.value)]
pigcell_mtcounts.value <- pigcell_mtcounts.value[,colnames(pigcell_counts.value)]
pigcell_stagecol <- pigcell_stagecol[colnames(pigcell_counts.value),]

colnames(pigcell_counts.value) <- pigcell_id[colnames(pigcell_counts.value),1]
colnames(pigcell_tpms.value) <- pigcell_id[colnames(pigcell_tpms.value),1]
colnames(pigcell_mtcounts.value) <- pigcell_id[colnames(pigcell_mtcounts.value),1]
names(pigcell_stagecol) <- pigcell_id[names(pigcell_stagecol),1]


#################################################################################################################################################
#                                                           separate cells by stage                                                             #
#################################################################################################################################################

pigcell_counts.value.e5_id <- names(pigcell_stagecol[pigcell_stagecol=="tomato"])
pigcell_counts.value.e6_id <- names(pigcell_stagecol[pigcell_stagecol=="skyblue"])
pigcell_counts.value.e8_id <- names(pigcell_stagecol[pigcell_stagecol=="yellowgreen"])
pigcell_counts.value.e11_id <- names(pigcell_stagecol[pigcell_stagecol=="orange"])
pigcell_counts.value.e14_id <- names(pigcell_stagecol[pigcell_stagecol=="purple"])
pigcell_counts.value.e31_id <- names(pigcell_stagecol[pigcell_stagecol=="blue"])

pigcell_counts.value.e5_11 <- pigcell_counts.value[,c(pigcell_counts.value.e5_id,pigcell_counts.value.e6_id,pigcell_counts.value.e8_id,pigcell_counts.value.e11_id)]
pigcell_counts.value.e14_31 <- pigcell_counts.value[,c(pigcell_counts.value.e14_id,pigcell_counts.value.e31_id)]


#################################################################################################################################################
#                                                           specific gene list                                                                  #
#################################################################################################################################################

genelist.e5_11 <- read.table("genelist_e5e6e8e11.txt", header = TRUE, row.names = 1)
genelist.e5_11 <- as.matrix(genelist.e5_11)

genelist.e14_31 <- read.table("genelist_e14e31.txt", header = TRUE, row.names = 1)
genelist.e14_31 <- as.matrix(genelist.e14_31)

genelist.e5_11_plur <- read.table("genelist_e5e6e8e11_pluripotency.txt", header = TRUE, row.names = 1)
genelist.e5_11_plur <- as.matrix(genelist.e5_11_plur)


#################################################################################################################################################
#                                                       filtered by total read counts                                                           #
#################################################################################################################################################

pigcell_counts.e5_11.countfilter <- pigcell_counts.value.e5_11[,colSums(pigcell_counts.value.e5_11)>1000000]
pigcell_counts.e14_31.countfilter <- pigcell_counts.value.e14_31[,colSums(pigcell_counts.value.e14_31)>1000000]


#################################################################################################################################################
#                                                   filter noise cells by mt expression ratio                                                   #
#################################################################################################################################################

pigcell_mtratio <- colSums(pigcell_mtcounts.value)/colSums(pigcell_counts.value[,colnames(pigcell_mtcounts.value)])

pigcell_mtratio.e5_11 <- pigcell_mtratio[colnames(pigcell_counts.value.e5_11)]
pigcell_mtratio.e5_11.countfilter <- pigcell_mtratio.e5_11[colnames(pigcell_counts.e5_11.countfilter)]
pigcell_counts.e5_11.mtfilter <- pigcell_counts.e5_11.countfilter[,names(pigcell_mtratio.e5_11.countfilter[pigcell_mtratio.e5_11.countfilter<0.5])]

pigcell_mtratio.e14_31 <- pigcell_mtratio[colnames(pigcell_counts.value.e14_31)]
pigcell_mtratio.e14_31.countfilter <- pigcell_mtratio.e14_31[colnames(pigcell_counts.e14_31.countfilter)]
pigcell_counts.e14_31.mtfilter <- pigcell_counts.e14_31.countfilter[,names(pigcell_mtratio.e14_31.countfilter[pigcell_mtratio.e14_31.countfilter<0.5])]


#################################################################################################################################################
#                                                               normalisation                                                                   #
#################################################################################################################################################

########
#### filter low expressed genes
########

####
# e5 e6 e8 e11
####

pigcell_sce.e5_11 <- newSCESet(countData=pigcell_counts.e5_11.mtfilter)
pigcell_sce.e5_11 <- calculateQCMetrics(pigcell_sce.e5_11)

# mean, calcAverage can caculate and normalize mean counts by library size (total counts)
pigcell_sce.e5_11.meancounts <- calcAverage(pigcell_sce.e5_11)
pigcell_sce.e5_11.meancounst_a1 <- pigcell_sce.e5_11.meancounts >= 0
sum(pigcell_sce.e5_11.meancounst_a1)

# gene expressed at least in 3 cells
pigcell_sce.e5_11.numcells <- nexprs(pigcell_sce.e5_11, byrow=TRUE)
pigcell_sce.e5_11.numcells_a10 <- pigcell_sce.e5_11.numcells >= 3
sum(pigcell_sce.e5_11.numcells_a10)

# remove genes with mean counts < 0, and remove genes expressed < 3 cells
pigcell_sce.e5_11.meancounst_a1_matrix <- pigcell_counts.e5_11.mtfilter[names(pigcell_sce.e5_11.meancounst_a1[pigcell_sce.e5_11.meancounst_a1]),]
pigcell_sce.e5_11.meancounst_a1_matrix_sce <- newSCESet(countData=pigcell_sce.e5_11.meancounst_a1_matrix)
pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells <- nexprs(pigcell_sce.e5_11.meancounst_a1_matrix_sce, byrow=TRUE)
pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells_a10 <- pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells >= 3
pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells_a10_matrix <- pigcell_sce.e5_11.meancounst_a1_matrix[names(pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells_a10[pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells_a10]),]

pigcell_sce.e5_11.genefilter <- pigcell_sce.e5_11.meancounst_a1_matrix_sce.numcells_a10_matrix
pigcell_stagecol.e5_11.genefilter <- pigcell_stagecol[colnames(pigcell_counts.e5_11.mtfilter)]

########
#### size factor (srcan) normalization
########

####
# e5 e6 e8 e11
####

pigcell.e5_11.genefilter_sce <- newSCESet(countData=pigcell_sce.e5_11.genefilter)
pigcell.e5_11.genefilter_sce <- calculateQCMetrics(pigcell.e5_11.genefilter_sce)
pigcell.e5_11.genefilter_sce <- computeSumFactors(pigcell.e5_11.genefilter_sce, sizes=c(1:10)*20)
summary(sizeFactors(pigcell.e5_11.genefilter_sce))
plot(sizeFactors(pigcell.e5_11.genefilter_sce), pigcell.e5_11.genefilter_sce$total_counts/1e6, log = "xy", ylab = "Library size (millions)", xlab = "Size factor", pch = 16, col = "brown1")
abline(lm(log10(pigcell.e5_11.genefilter_sce$total_counts/1e6) ~ log10(sizeFactors(pigcell.e5_11.genefilter_sce))), col = "blue", lty = 2, lwd = 2)

# notice: exprs return log2(counts-per-million) using the cpm function from edgeR, prior count = 1
pigcell.e5_11.genefilter_sce <- normalize(pigcell.e5_11.genefilter_sce)
pigcell.e5_11.genefilter_sce.norm_exprs <- exprs(pigcell.e5_11.genefilter_sce)

#################################################################################################################################################
#                                                                   Cluster Trees                                                               #
#################################################################################################################################################

########
#### Ward Hierarchical Clustering -- cells
########

####
# e5 e6 e8 e11
####

pigcell.e5_11.genefilter_sce.dist <- dist(t(pigcell.e5_11.genefilter_sce.norm_exprs), method = "euclidean")
pigcell.e5_11.genefilter_sce.hc <- hclust(pigcell.e5_11.genefilter_sce.dist, method="ward.D2")

pigcell.e5_11.genefilter_sce.dend <- as.dendrogram(pigcell.e5_11.genefilter_sce.hc)

pigcell.e5_11.genefilter_sce.dend <- set(pigcell.e5_11.genefilter_sce.dend, "by_labels_branches_col", value = names(pigcell_stagecol.e5_11.genefilter[pigcell_stagecol.e5_11.genefilter=="tomato"]), TF_value = "tomato", type = "any")
pigcell.e5_11.genefilter_sce.dend <- set(pigcell.e5_11.genefilter_sce.dend, "by_labels_branches_col", value = names(pigcell_stagecol.e5_11.genefilter[pigcell_stagecol.e5_11.genefilter=="skyblue"]), TF_value = "skyblue", type = "any")
pigcell.e5_11.genefilter_sce.dend <- set(pigcell.e5_11.genefilter_sce.dend, "by_labels_branches_col", value = names(pigcell_stagecol.e5_11.genefilter[pigcell_stagecol.e5_11.genefilter=="yellowgreen"]), TF_value = "yellowgreen", type = "any")
pigcell.e5_11.genefilter_sce.dend <- set(pigcell.e5_11.genefilter_sce.dend, "by_labels_branches_col", value = names(pigcell_stagecol.e5_11.genefilter[pigcell_stagecol.e5_11.genefilter=="orange"]), TF_value = "orange", type = "any")

plot(pigcell.e5_11.genefilter_sce.dend, leaf = "none", horiz = TRUE)
rect.dendrogram(pigcell.e5_11.genefilter_sce.dend, h = 1000, horiz = TRUE, border = "darkgrey", lty = 2, lwd = 1)

par(mar=c(4,4,1,10), mfrow=c(2,2), cex = 0.4)
plot(cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[1]], horiz = T)
plot(cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[2]], horiz = T)
plot(cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[3]], horiz = T)
plot(cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[4]], horiz = T)

# tree
pigcell_counts.e5_11.cut <- cutree(pigcell.e5_11.genefilter_sce.dend, h = 650)
pigcell_counts.e5_11.h650_tree1_e5 <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==1]
pigcell_counts.e5_11.h650_tree2_e6icm <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==2]
pigcell_counts.e5_11.h650_tree3_e6te <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==3]
pigcell_counts.e5_11.h650_tree4_e8hypo <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==4]
pigcell_counts.e5_11.h650_tree5_e8epi <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==5]
pigcell_counts.e5_11.h650_tree6_e11hypo <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==6]
pigcell_counts.e5_11.h650_tree7_e11epi <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==7]

pigcell_counts.e5_11.cut <- cutree(pigcell.e5_11.genefilter_sce.dend, h = 470)
pigcell_counts.e5_11.h470_tree8_e8epi_g1 <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==8]
pigcell_counts.e5_11.h470_tree10_e8epi_g2 <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==10]
pigcell_counts.e5_11.h470_tree9_e8hypo_g3 <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==9]
pigcell_counts.e5_11.h470_tree7_e8hypo_g4 <- pigcell_counts.e5_11.cut[pigcell_counts.e5_11.cut==7]

# dendrogram
pigcell_counts.e5_11.h650_dend1_e8epi <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[1]]
pigcell_counts.e5_11.h650_dend2_e11epi <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[2]]
pigcell_counts.e5_11.h650_dend3_e8hypo <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[3]]
pigcell_counts.e5_11.h650_dend4_e11hypo <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[4]]
pigcell_counts.e5_11.h650_dend5_e5 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[5]]
pigcell_counts.e5_11.h650_dend6_e6icm <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[6]]
pigcell_counts.e5_11.h650_dend7_e6te <- cut(pigcell.e5_11.genefilter_sce.dend, h = 650)$lower[[7]]

pigcell_counts.e5_11.h800_dend1_e8epi_e11epi <- cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[1]]
pigcell_counts.e5_11.h800_dend2_e8hypo_e11hypo <- cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[2]]
pigcell_counts.e5_11.h800_dend3_e5 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[3]]
pigcell_counts.e5_11.h800_dend4_e6 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 800)$lower[[4]]

pigcell_counts.e5_11.h470_dend1_e8epi_g1 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 470)$lower[[1]]
pigcell_counts.e5_11.h470_dend2_e8epi_g2 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 470)$lower[[2]]
pigcell_counts.e5_11.h470_dend5_e8hypo_g3 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 470)$lower[[5]]
pigcell_counts.e5_11.h470_dend6_e8hypo_g4 <- cut(pigcell.e5_11.genefilter_sce.dend, h = 470)$lower[[6]]

# set dendrogram color
pigcell_counts.e5_11.h800_dend1_e8epi_e11epi <- set(pigcell_counts.e5_11.h800_dend1_e8epi_e11epi, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend1_e8epi), TF_value = "lightgreen", type = "any")
pigcell_counts.e5_11.h800_dend1_e8epi_e11epi <- set(pigcell_counts.e5_11.h800_dend1_e8epi_e11epi, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend2_e11epi), TF_value = "coral", type = "any")

pigcell_counts.e5_11.h800_dend2_e8hypo_e11hypo <- set(pigcell_counts.e5_11.h800_dend2_e8hypo_e11hypo, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend3_e8hypo), TF_value = "yellowgreen", type = "any")
pigcell_counts.e5_11.h800_dend2_e8hypo_e11hypo <- set(pigcell_counts.e5_11.h800_dend2_e8hypo_e11hypo, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend4_e11hypo), TF_value = "orange", type = "any")

pigcell_counts.e5_11.h800_dend4_e6 <- set(pigcell_counts.e5_11.h800_dend4_e6, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend6_e6icm), TF_value = "royalblue1", type = "any")
pigcell_counts.e5_11.h800_dend4_e6 <- set(pigcell_counts.e5_11.h800_dend4_e6, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend7_e6te), TF_value = "skyblue", type = "any")

pigcell_counts.e5_11.h650_dend1_e8epi <- set(pigcell_counts.e5_11.h650_dend1_e8epi, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend1_e8epi), TF_value = "lightgreen", type = "any")
pigcell_counts.e5_11.h650_dend3_e8hypo <- set(pigcell_counts.e5_11.h650_dend3_e8hypo, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend3_e8hypo), TF_value = "yellowgreen", type = "any")

pigcell_counts.e5_11.h650_dend2_e11epi <- set(pigcell_counts.e5_11.h650_dend2_e11epi, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend2_e11epi), TF_value = "coral", type = "any")
pigcell_counts.e5_11.h650_dend4_e11hypo <- set(pigcell_counts.e5_11.h650_dend4_e11hypo, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h650_dend4_e11hypo), TF_value = "gold", type = "any")

pigcell_counts.e5_11.h470_dend1_e8epi_g1 <- set(pigcell_counts.e5_11.h470_dend1_e8epi_g1, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h470_dend1_e8epi_g1), TF_value = "lightgreen", type = "any")
pigcell_counts.e5_11.h470_dend2_e8epi_g2 <- set(pigcell_counts.e5_11.h470_dend2_e8epi_g2, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h470_dend2_e8epi_g2), TF_value = "lightgreen", type = "any")
pigcell_counts.e5_11.h470_dend5_e8hypo_g3 <- set(pigcell_counts.e5_11.h470_dend5_e8hypo_g3, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h470_dend5_e8hypo_g3), TF_value = "yellowgreen", type = "any")
pigcell_counts.e5_11.h470_dend6_e8hypo_g4 <- set(pigcell_counts.e5_11.h470_dend6_e8hypo_g4, "by_labels_branches_col", value = labels(pigcell_counts.e5_11.h470_dend6_e8hypo_g4), TF_value = "yellowgreen", type = "any")

###########################################################################################################################################
#                                                                   PCA                                                                   #
###########################################################################################################################################

########
#### PCA & Loadings
########

e5_e6icm_e8epi_e11epi.cells <- c(labels(pigcell_counts.e5_11.h650_dend5_e5), labels(pigcell_counts.e5_11.h650_dend6_e6icm), labels(pigcell_counts.e5_11.h650_dend1_e8epi), labels(pigcell_counts.e5_11.h650_dend2_e11epi))
e5_e6icm_e8epi_e11epi.tpm <- pigcell_tpms.value[,e5_e6icm_e8epi_e11epi.cells]

e5_e6icm_e8epi_e11epi.tpm.mean <- rowMeans(e5_e6icm_e8epi_e11epi.tpm)
e5_e6icm_e8epi_e11epi.tpm.mean.a0 <- e5_e6icm_e8epi_e11epi.tpm[names(e5_e6icm_e8epi_e11epi.tpm.mean[e5_e6icm_e8epi_e11epi.tpm.mean>0]),]
e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10 <- log10(e5_e6icm_e8epi_e11epi.tpm.mean.a0+1)

e5_e6icm_e8epi_e11epi.col <- pigcell_stagecol.e5_11.genefilter[e5_e6icm_e8epi_e11epi.cells]
e5_e6icm_e8epi_e11epi.col[names(pigcell_counts.e5_11.h650_tree1_e5)] <- "red"
e5_e6icm_e8epi_e11epi.col[names(pigcell_counts.e5_11.h650_tree2_e6icm)] <- "royalblue1"
e5_e6icm_e8epi_e11epi.col[names(pigcell_counts.e5_11.h650_tree5_e8epi)] <- "lightgreen"
e5_e6icm_e8epi_e11epi.col[names(pigcell_counts.e5_11.h650_tree7_e11epi)] <- "coral"

e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca <- prcomp(t(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10), center = TRUE, scale = TRUE)

par(mfrow = c(2,2), mar = c(5,5,1,1))
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,1], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,2], xlab = "PC1 (6.40%)", ylab = "PC2 (2.44%)", pch = 21, bg = e5_e6icm_e8epi_e11epi.col, col = "black", lwd = 0.3, frame = F, xlim = c(-60,60), ylim = c(-100,100), cex.lab = 2, cex = 2)
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,1], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,3], xlab = "PC1 (6.40%)", ylab = "PC3 (2.33%)", pch = 21, bg = e5_e6icm_e8epi_e11epi.col, col = "black", lwd = 0.3, frame = F, xlim = c(-60,60), ylim = c(-60,60), cex.lab = 2, cex = 2)
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,2], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,3], xlab = "PC2 (2.44%)", ylab = "PC3 (2.33%)", pch = 21, bg = e5_e6icm_e8epi_e11epi.col, col = "black", lwd = 0.3, frame = F, xlim = c(-100,100), ylim = c(-60,60), cex.lab = 2, cex = 2)

scatterplot3d(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,3], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,2], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$x[,1], xlab = "PC3 (2.33%)", ylab = "PC2 (2.44%)", zlab = "PC1 (6.40%)", box = F, pch = 21, bg = e5_e6icm_e8epi_e11epi.col, lwd = 0.1, angle = 15, cex.symbols = 1.5, cex.lab = 2, cex.axis = 1, xlim = c(-60,60), ylim = c(-60,60), zlim = c(-60,60))

e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings <- scale(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$rotation[,1:3])
e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col <- scale(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca$rotation[,1:3])

pc1.sd <- sd(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,1])
pc2.sd <- sd(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,2])
pc3.sd <- sd(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,3])

for(i in 1:dim(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings)[1])
{
    #### two PCs
    
    if((e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,1])^2 + (e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,2])^2 > (3*pc1.sd)^2)
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,1] <- "brown1"
    }
    else
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,1] <- "lightgrey"
    }
    
    if((e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,1])^2 + (e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,3])^2 > (3*pc1.sd)^2)
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,2] <- "brown1"
    }
    else
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,2] <- "lightgrey"
    }
    
    if((e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,2])^2 + (e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[i,3])^2 > (3*pc1.sd)^2)
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,3] <- "brown1"
    }
    else
    {
        e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[i,3] <- "lightgrey"
    }
}

par(mfrow = c(2,2), mar = c(5,5,1,1))
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,1], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,2], xlab = "z score of PC1", ylab = "z score of PC2", pch = 21, bg = e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,1], col = "darkgrey", lwd = 0.1, xlim = c(-6,6), ylim = c(-6,6), cex = 1, cex.lab = 2, frame = F)
text(-2, 6, "rPC12 > SD3: 436 genes", cex = 1.5)
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,1], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,3], xlab = "z score of PC1", ylab = "z score of PC3", pch = 21, bg = e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,2], col = "darkgrey", lwd = 0.1, xlim = c(-6,6), ylim = c(-6,6), cex = 1, cex.lab = 2, frame = F)
text(-2, 6, "rPC13 > SD3: 564 genes", cex = 1.5)
plot(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,2], e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings[,3], xlab = "z score of PC2", ylab = "z score of PC3", pch = 21, bg = e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,3], col = "darkgrey", lwd = 0.1, xlim = c(-6,6), ylim = c(-6,6), cex = 1, cex.lab = 2, frame = F)
text(-2, 6, "rPC23 > SD3: 393 genes", cex = 1.5)

pc12.genes <- rownames(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,1]=="brown1",])
e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pc12sig <- e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10[pc12.genes,]

pc13.genes <- rownames(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,2]=="brown1",])
e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pc13sig <- e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10[pc13.genes,]

pc23.genes <- rownames(e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pca.loadings.col[,3]=="brown1",])
e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10.pc23sig <- e5_e6icm_e8epi_e11epi.tpm.mean.a0.log10[pc23.genes,]

#################################################################################################################################################
#                                                                Pig Monkey Mouse                                                               #
#################################################################################################################################################

########
#### PCA & tSNE
########

pig_monkey_mouse.rpms <- read.table("ensembl.homologies.Compara.87.pig_monkey_mouse.ncbi.unique.rpms.txt", header = T, row.names = 2, sep = '\t')
pig_monkey_mouse.rpms.value <- pig_monkey_mouse.rpms[,4:dim(pig_monkey_mouse.rpms)[2]]
colnames(pig_monkey_mouse.rpms.value)[1:393] <- pigcell_id[colnames(pig_monkey_mouse.rpms.value)[1:393],]

pig.cond1.cells <- c(labels(pigcell_counts.e5_11.h650_dend6_e6icm),labels(pigcell_counts.e5_11.h650_dend1_e8epi),labels(pigcell_counts.e5_11.h650_dend2_e11epi))
monkey.cond1.cells <- read.table("macFas5.selectedcell.cond1.txt", row.names = 1)
monkey.cond1.cells <- as.matrix(monkey.cond1.cells)
mouse.cond1.cells <- read.table("mm10.selectedcell.cond1.txt", row.names = 1)
mouse.cond1.cells <- as.matrix(mouse.cond1.cells)

pig_monkey_mouse.rpms.cond1 <- pig_monkey_mouse.rpms.value[,c(pig.cond1.cells,rownames(monkey.cond1.cells),rownames(mouse.cond1.cells))]
pig_monkey_mouse.rpms.cond1.mean <- rowMeans(pig_monkey_mouse.rpms.cond1)
pig_monkey_mouse.rpms.cond1.a0 <- pig_monkey_mouse.rpms.cond1[names(pig_monkey_mouse.rpms.cond1.mean[pig_monkey_mouse.rpms.cond1.mean!=0]),]

pig_monkey_mouse.col <- c(rep("royalblue1", length(labels(pigcell_counts.e5_11.h650_dend6_e6icm))),
                          rep("lightgreen", length(labels(pigcell_counts.e5_11.h650_dend1_e8epi))),
                          rep("coral", length(labels(pigcell_counts.e5_11.h650_dend2_e11epi))),
                          rep("slateblue", sum(monkey.cond1.cells=="ICM")),
                          rep("yellowgreen", sum(monkey.cond1.cells=="Pre-EPI")),
                          rep("tomato", sum(monkey.cond1.cells=="PostE-EPI")),
                          rep("skyblue", sum(mouse.cond1.cells=="E3.5ICM")),
                          rep("seagreen1", sum(mouse.cond1.cells=="E4.5EPI")),
                          rep("orange", sum(mouse.cond1.cells=="E5.5EPI")))
names(pig_monkey_mouse.col) <- colnames(pig_monkey_mouse.rpms.cond1)

pig_monkey_mouse.pch <- c(rep(21, length(labels(pigcell_counts.e5_11.h650_dend6_e6icm))),
                          rep(21, length(labels(pigcell_counts.e5_11.h650_dend1_e8epi))),
                          rep(21, length(labels(pigcell_counts.e5_11.h650_dend2_e11epi))),
                          rep(23, sum(monkey.cond1.cells=="ICM")),
                          rep(23, sum(monkey.cond1.cells=="Pre-EPI")),
                          rep(23, sum(monkey.cond1.cells=="PostE-EPI")),
                          rep(24, sum(mouse.cond1.cells=="E3.5ICM")),
                          rep(24, sum(mouse.cond1.cells=="E4.5EPI")),
                          rep(24, sum(mouse.cond1.cells=="E5.5EPI")))
names(pig_monkey_mouse.pch) <- colnames(pig_monkey_mouse.rpms.cond1)

####
# PCA
####

pig_monkey_mouse.pca <- prcomp(t(log2(pig_monkey_mouse.rpms.cond1.a0+1)), center = TRUE, scale = TRUE)
summary(pig_monkey_mouse.pca)

plot3d(pig_monkey_mouse.pca$x[,1], pig_monkey_mouse.pca$x[,2], pig_monkey_mouse.pca$x[,3], type = 'p', size = 10, col = pig_monkey_mouse.col)

par(mfrow = c(2,3), mar = c(4,5,1,1))

plot(pig_monkey_mouse.pca$x[,1], pig_monkey_mouse.pca$x[,2], xlab = "PC1 (16.73%)", ylab = "PC2 (10.35%)", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F, xlim = c(-100,100), ylim = c(-100,100))
plot(pig_monkey_mouse.pca$x[,1], pig_monkey_mouse.pca$x[,3], xlab = "PC1 (16.73%)", ylab = "PC3 (2.63%)", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F, xlim = c(-100,100), ylim = c(-60,60))
plot(pig_monkey_mouse.pca$x[,2], pig_monkey_mouse.pca$x[,3], xlab = "PC2 (10.35%)", ylab = "PC3 (2.63%)", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F, xlim = c(-100,100), ylim = c(-60,60))

par(mar = c(0,0,0,0))
plot.new()
legend("center", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, cex = 1.5, pch = c(rep(19,3),rep(18,3),rep(17,3)), bty = "n")

scatterplot3d(pig_monkey_mouse.pca$x[,1], pig_monkey_mouse.pca$x[,2], pig_monkey_mouse.pca$x[,3], xlab = "PC1 (16.73%)", ylab = "PC2 (10.35%)", zlab = "PC3 (2.63%)", box = F, pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, lwd = 0.2, angle = 10, cex.symbols = 1.5, cex.lab = 1.5, cex.axis = 1)

plot3d(pig_monkey_mouse.pca $x[,1], pig_monkey_mouse.pca $x[,3], pig_monkey_mouse.pca $x[,2], type = 'p', size = 10, col = pig_monkey_mouse.col, box = FALSE, xlab = "PC1 (16.73%)", ylab = "PC3 (2.63%)", zlab = "PC2 (10.35%)", expand = 1.2)
text3d(c(-55,-20,60), c(15,15,0), c(-40,60,-10), c("Pig","Mouse","Monkey"), col = "red")
legend3d("topright", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, pch = 16, cex = 1.5, bty = "n")

####
# tSNE
####

library(Rtsne)

pig_monkey_mouse.tsne <- Rtsne(t(log2(pig_monkey_mouse.rpms.cond1.a0+1)), dim = 3, perplexity = 20, max_iter = 2000)

plot3d(pig_monkey_mouse.tsne$Y, type = 'p', size = 10, col = pig_monkey_mouse.col)

par(mfrow = c(2,3), mar = c(5,5,1,1))
plot(pig_monkey_mouse.tsne$Y[,1], pig_monkey_mouse.tsne$Y[,2], xlab = "tSNE1", ylab = "tSNE2", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F)
plot(pig_monkey_mouse.tsne$Y[,1], pig_monkey_mouse.tsne$Y[,3], xlab = "tSNE1", ylab = "tSNE3", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F)
plot(pig_monkey_mouse.tsne$Y[,2], pig_monkey_mouse.tsne$Y[,3], xlab = "tSNE2", ylab = "tSNE3", pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, col = "black", lwd = 0.2, cex = 2, cex.lab = 2, frame = F)

par(mar = c(0,0,0,0))
plot.new()
legend("center", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, cex = 1.5, pch = c(rep(19,3),rep(18,3),rep(17,3)), bty = "n")

scatterplot3d(pig_monkey_mouse.tsne$Y[,2], pig_monkey_mouse.tsne$Y[,3], pig_monkey_mouse.tsne$Y[,1], xlab = "tSNE2", ylab = "tSNE3", zlab = "tSNE1", box = F, pch = pig_monkey_mouse.pch, bg = pig_monkey_mouse.col, lwd = 0.2, angle = 10, cex.symbols = 1.5, cex.lab = 1.5, cex.axis = 1)

########
#### HeatMap
########

pig_monkey_mouse.pca.loadings <- scale(pig_monkey_mouse.pca$rotation[,1:3])
pig_monkey_mouse.pca.loadings.col <- scale(pig_monkey_mouse.pca$rotation[,1:6])

pc1.sd <- sd(pig_monkey_mouse.pca.loadings[,1])
pc2.sd <- sd(pig_monkey_mouse.pca.loadings[,2])
pc3.sd <- sd(pig_monkey_mouse.pca.loadings[,3])

for(i in 1:dim(pig_monkey_mouse.pca.loadings)[1])
{
    #### single PC
    
    if(abs(pig_monkey_mouse.pca.loadings[i,1]) > 2*pc1.sd)
    {
        pig_monkey_mouse.pca.loadings.col[i,1] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,1] <- "lightgrey"
    }
    
    if(abs(pig_monkey_mouse.pca.loadings[i,2]) > 2*pc1.sd)
    {
        pig_monkey_mouse.pca.loadings.col[i,2] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,2] <- "lightgrey"
    }
    
    if(abs(pig_monkey_mouse.pca.loadings[i,3]) > 2*pc1.sd)
    {
        pig_monkey_mouse.pca.loadings.col[i,3] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,3] <- "lightgrey"
    }
    
    #### double PCs
    
    if((pig_monkey_mouse.pca.loadings[i,1])^2 + (pig_monkey_mouse.pca.loadings[i,2])^2 > (3*pc1.sd)^2)
    {
        pig_monkey_mouse.pca.loadings.col[i,4] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,4] <- "lightgrey"
    }
    
    if((pig_monkey_mouse.pca.loadings[i,1])^2 + (pig_monkey_mouse.pca.loadings[i,3])^2 > (3*pc1.sd)^2)
    {
        pig_monkey_mouse.pca.loadings.col[i,5] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,5] <- "lightgrey"
    }
    
    if((pig_monkey_mouse.pca.loadings[i,2])^2 + (pig_monkey_mouse.pca.loadings[i,3])^2 > (3*pc1.sd)^2)
    {
        pig_monkey_mouse.pca.loadings.col[i,6] <- "brown1"
    }
    else
    {
        pig_monkey_mouse.pca.loadings.col[i,6] <- "lightgrey"
    }
}

pc1.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,1]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc1sig <- pig_monkey_mouse.rpms.cond1.a0[pc1.genes,]

pc2.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,2]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc2sig <- pig_monkey_mouse.rpms.cond1.a0[pc2.genes,]

pc3.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,3]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc3sig <- pig_monkey_mouse.rpms.cond1.a0[pc3.genes,]

pc12.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,4]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc12sig <- pig_monkey_mouse.rpms.cond1.a0[pc12.genes,]

pc13.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,5]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc13sig <- pig_monkey_mouse.rpms.cond1.a0[pc13.genes,]

pc23.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,6]=="brown1",])
pig_monkey_mouse.rpms.cond1.a0.pc23sig <- pig_monkey_mouse.rpms.cond1.a0[pc23.genes,]

####
# pig_monkey_mouse pca sig PC1
####

rownames(pig_monkey_mouse.rpms.cond1.a0.pc1sig) <- rownames(pigcell.geneid[pigcell.geneid[,1]%in%rownames(pig_monkey_mouse.rpms.cond1.a0.pc1sig),])

lmat = rbind(c(5,0,0),c(0,4,1),c(6,3,2))
lhei = c(0.9,0.2,8)
lwid = c(3,2,6)

layout(lmat, widths = lwid, heights = lhei)

pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm <- heatmap.2.nolayout(as.matrix(log2(pig_monkey_mouse.rpms.cond1.a0.pc1sig+1)),
                                                               #hclustfun = function(x) hclust(x,method = 'ward.D2'),
                                                               col = "bluered",
                                                               density.info = "none",
                                                               trace = "none",
                                                               dendrogram = "row",
                                                               Colv = FALSE,
                                                               ColSideColors = pig_monkey_mouse.col,
                                                               labRow = F, labCol = F,
                                                               key.xlab = "log2(rpm+1)",
                                                               key.title = "",
                                                               key.par = list(cex.lab = 1),
                                                               margin = c(1,1))

par(mar = c(0,0,0,0))
plot.new()
legend("top", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, cex = 1.5, pch = c(rep(19,3),rep(18,3),rep(17,3)), bty = "n")

pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[2]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend3 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[3]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend4 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[4]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend5 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[5]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend6 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[6]]),]
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend7 <- pig_monkey_mouse.rpms.cond1.a0.pc1sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc1sig.hm$rowDendrogram, h = 100)$lower[[7]]),]

pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend2")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend3.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend3, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend3")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend4.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend4, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend4")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend5.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend5, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend5")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend6.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend6, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend6")
pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend7.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend7, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc1sig", "pig_monkey_mouse.rpms.cond1.a0.pc1sig.h100_dend7")

####
# pig_monkey_mouse pca sig PC2
####

rownames(pig_monkey_mouse.rpms.cond1.a0.pc2sig) <- rownames(pigcell.geneid[pigcell.geneid[,1]%in%rownames(pig_monkey_mouse.rpms.cond1.a0.pc2sig),])

lmat = rbind(c(5,0,0),c(0,4,1),c(6,3,2))
lhei = c(0.9,0.2,8)
lwid = c(3,2,6)

layout(lmat, widths = lwid, heights = lhei)

pig_monkey_mouse.rpms.cond1.a0.pc2sig.hm <- heatmap.2.nolayout(as.matrix(log2(pig_monkey_mouse.rpms.cond1.a0.pc2sig+1)),
                                                               #hclustfun = function(x) hclust(x,method = 'ward.D2'),
                                                               col = "bluered",
                                                               density.info = "none",
                                                               trace = "none",
                                                               dendrogram = "row",
                                                               Colv = FALSE,
                                                               ColSideColors = pig_monkey_mouse.col,
                                                               labRow = F, labCol = F,
                                                               key.xlab = "log2(rpm+1)",
                                                               key.title = "",
                                                               key.par = list(cex.lab = 1),
                                                               margin = c(1,1))

par(mar = c(0,0,0,0))
plot.new()
legend("top", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, cex = 1.5, pch = c(rep(19,3),rep(18,3),rep(17,3)), bty = "n")

pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc2sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.hm$rowDendrogram, h = 150)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc2sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.hm$rowDendrogram, h = 150)$lower[[2]]),]
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend3 <- pig_monkey_mouse.rpms.cond1.a0.pc2sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.hm$rowDendrogram, h = 150)$lower[[3]]),]

pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.dend <- cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.hm$rowDendrogram, h = 150)$lower[[2]]
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc2sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.dend, h = 100)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc2sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.dend, h = 100)$lower[[2]]),]

pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2.h100_dend2")
pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend3.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend3, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend3")

pig_monkey_mouse.rpms.cond1.a0.pc2sig.up <- rbind(pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend2,pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend3)
pig_monkey_mouse.rpms.cond1.a0.pc2sig.down <- pig_monkey_mouse.rpms.cond1.a0.pc2sig.h150_dend1

pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.up, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.up")
pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc2sig.down, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc2sig", "pig_monkey_mouse.rpms.cond1.a0.pc2sig.down")

pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.bp_go <- pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.go[pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.go[,7]=="BP"&(!is.na(pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.go[,7])),]
pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.bp_go <- pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.go[pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.go[,7]=="BP"&(!is.na(pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.go[,7])),]

par(mfrow=c(1,2))

par(mar=c(5,40,2,1.1))
barplot(rev(-log10(pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.bp_go[1:30,2])), horiz = TRUE, names.arg = rev(pig_monkey_mouse.rpms.cond1.a0.pc2sig.up.bp_go[1:30,6]), las = 2, xlim = rev(c(0,8)), col = rgb(255, 30, 0, alpha=160, maxColorValue=255), border = NA, xlab = "-log10(p-value)", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5, main = "Up-regualted", col.main = "brown1")

par(mar=c(5,1.1,2,40))
m <- barplot(rev(-log10(pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.bp_go[1:30,2])), horiz = TRUE, names.arg = NA, las = 2, xlim = c(0,6), col = rgb(0, 0, 255, alpha=160, maxColorValue=255), border = NA, xlab = "-log10(p-value)", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5, main = "Down-regulated", col.main = "slateblue")
axis(4, m[c(1:length(pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.bp_go[1:30,2]))], labels = rev(pig_monkey_mouse.rpms.cond1.a0.pc2sig.down.bp_go[1:30,6]), las = 2, tick = FALSE, cex.axis = 1.5)

####
# pig_monkey_mouse pca sig PC3
####

rownames(pig_monkey_mouse.rpms.cond1.a0.pc3sig) <- rownames(pigcell.geneid[pigcell.geneid[,1]%in%rownames(pig_monkey_mouse.rpms.cond1.a0.pc3sig),])

lmat = rbind(c(5,0,0),c(0,4,1),c(6,3,2))
lhei = c(0.9,0.2,8)
lwid = c(3,2,6)

layout(lmat, widths = lwid, heights = lhei)

pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm <- heatmap.2.nolayout(as.matrix(log2(pig_monkey_mouse.rpms.cond1.a0.pc3sig+1)),
                                                               hclustfun = function(x) hclust(x,method = 'ward.D2'),
                                                               col = "bluered",
                                                               density.info = "none",
                                                               trace = "none",
                                                               dendrogram = "row",
                                                               Colv = FALSE,
                                                               ColSideColors = pig_monkey_mouse.col,
                                                               labRow = F, labCol = F,
                                                               key.xlab = "log2(rpm+1)",
                                                               key.title = "",
                                                               key.par = list(cex.lab = 1),
                                                               margin = c(1,1))

par(mar = c(0,0,0,0))
plot.new()
legend("top", inset = 0, legend = c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EP"), col = c("royalblue1","lightgreen","coral","slateblue","yellowgreen","tomato","skyblue","seagreen1","orange"), x.intersp = 1, cex = 1.5, pch = c(rep(19,3),rep(18,3),rep(17,3)), bty = "n")

pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[2]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend3 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[3]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[4]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[5]]),]

pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.dend <- cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[1]]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.dend, h = 300)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.dend, h = 300)$lower[[2]]),]

pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.dend <- cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[4]]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.dend, h = 200)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.dend, h = 200)$lower[[2]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend3 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.dend, h = 200)$lower[[3]]),]

pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.dend <- cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.hm$rowDendrogram, h = 400)$lower[[5]]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend1 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.dend, h = 300)$lower[[1]]),]
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend2 <- pig_monkey_mouse.rpms.cond1.a0.pc3sig[labels(cut(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.dend, h = 300)$lower[[2]]),]

pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend1.h300_dend2")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend2")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend3.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend3, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend3")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend2")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend3.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend3, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend4.h200_dend3")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend1.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend1, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend1")
pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend2.go <- run.goseq(pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend2, pig_genelen, pig_goterm, pigcell.geneid, "goseq_pc3sig", "pig_monkey_mouse.rpms.cond1.a0.pc3sig.h400_dend5.h300_dend2")

####
# pig_monkey_mouse pca sig PC23 and non-sig PC1
####

pc1.genes <- rownames(pig_monkey_mouse.pca.loadings.col[pig_monkey_mouse.pca.loadings.col[,1]=="lightgrey",])
pig_monkey_mouse.rpms.cond1.a0.pc1nonsig <- pig_monkey_mouse.rpms.cond1.a0[pc1.genes,]

pig_monkey_mouse.rpms.cond1.a0.pc23sig.pc1nonsig <- pig_monkey_mouse.rpms.cond1.a0.pc23sig[rownames(pig_monkey_mouse.rpms.cond1.a0.pc23sig)%in%rownames(pig_monkey_mouse.rpms.cond1.a0.pc1nonsig),]

rownames(pig_monkey_mouse.rpms.cond1.a0.pc23sig.pc1nonsig) <- rownames(pigcell.geneid[pigcell.geneid[,1]%in%rownames(pig_monkey_mouse.rpms.cond1.a0.pc23sig.pc1nonsig),])

lmat = rbind(c(5,0,0),c(0,4,1),c(6,3,2))
lhei = c(0.9,0.2,8)
lwid = c(3,2,6)

layout(lmat, widths = lwid, heights = lhei)

pig_monkey_mouse.rpms.cond1.a0.pc23sig.pc1nonsig.hm <- heatmap.2.nolayout(as.matrix(log2(pig_monkey_mouse.rpms.cond1.a0.pc23sig.pc1nonsig+1)),
                                                                          col = "bluered",
                                                                          density.info = "none",
                                                                          trace = "none",
                                                                          dendrogram = "row",
                                                                          Colv = FALSE,
                                                                          ColSideColors = pig_monkey_mouse.col,
                                                                          labRow = F, labCol = F,
                                                                          key.xlab = "log2(rpm+1)",
                                                                          key.title = "",
                                                                          key.par = list(cex.lab = 1),
                                                                          margin = c(1,1))

par(mar = c(0,0,0,0))
plot.new()
legend("top", inset = 0, legend = c("Pig_E5","Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Monkey_PostL-EPI","Mouse_E4.5EPI","Mouse_E5.5EPI","Mouse_E6.5EPI_T_lo","Mouse_E6.5EPI_T_hi"), col = c("red","royalblue1","lightgreen","coral","blue","forestgreen","yellowgreen","orange","pink","violet","skyblue","slateblue"), x.intersp = 1, cex = 1.5, pch = 15, bty = "n")

########
#### correlation coefficient
########

rowGeoMeans <- function(demat)
{
    outmat <- exp(rowMeans(log(demat+1)))
    
    return(outmat)
}

pig_monkey_mouse.rpms.cond1.a0.pig_e6icm <- pig_monkey_mouse.rpms.cond1.a0[,labels(pigcell_counts.e5_11.h650_dend6_e6icm)]
pig_monkey_mouse.rpms.cond1.a0.pig_e8epi <- pig_monkey_mouse.rpms.cond1.a0[,labels(pigcell_counts.e5_11.h650_dend1_e8epi)]
pig_monkey_mouse.rpms.cond1.a0.pig_e11epi <- pig_monkey_mouse.rpms.cond1.a0[,labels(pigcell_counts.e5_11.h650_dend2_e11epi)]
pig_monkey_mouse.rpms.cond1.a0.monkey_icm <- pig_monkey_mouse.rpms.cond1.a0[,names(monkey.cond1.cells[monkey.cond1.cells[,1]=="ICM",])]
pig_monkey_mouse.rpms.cond1.a0.monkey_pre_epi <- pig_monkey_mouse.rpms.cond1.a0[,names(monkey.cond1.cells[monkey.cond1.cells[,1]=="Pre-EPI",])]
pig_monkey_mouse.rpms.cond1.a0.monkey_poste_epi <- pig_monkey_mouse.rpms.cond1.a0[,names(monkey.cond1.cells[monkey.cond1.cells[,1]=="PostE-EPI",])]
pig_monkey_mouse.rpms.cond1.a0.mouse_e35icm <- pig_monkey_mouse.rpms.cond1.a0[,names(mouse.cond1.cells[mouse.cond1.cells[,1]=="E3.5ICM",])]
pig_monkey_mouse.rpms.cond1.a0.mouse_e45epi <- pig_monkey_mouse.rpms.cond1.a0[,names(mouse.cond1.cells[mouse.cond1.cells[,1]=="E4.5EPI",])]
pig_monkey_mouse.rpms.cond1.a0.mouse_e55epi <- pig_monkey_mouse.rpms.cond1.a0[,names(mouse.cond1.cells[mouse.cond1.cells[,1]=="E5.5EPI",])]

pig_monkey_mouse.rpms.cond1.a0.pig_e6icm.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.pig_e6icm)
pig_monkey_mouse.rpms.cond1.a0.pig_e8epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.pig_e8epi)
pig_monkey_mouse.rpms.cond1.a0.pig_e11epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.pig_e11epi)
pig_monkey_mouse.rpms.cond1.a0.monkey_icm.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.monkey_icm)
pig_monkey_mouse.rpms.cond1.a0.monkey_pre_epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.monkey_pre_epi)
pig_monkey_mouse.rpms.cond1.a0.monkey_poste_epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.monkey_poste_epi)
pig_monkey_mouse.rpms.cond1.a0.mouse_e35icm.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.mouse_e35icm)
pig_monkey_mouse.rpms.cond1.a0.mouse_e45epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.mouse_e45epi)
pig_monkey_mouse.rpms.cond1.a0.mouse_e55epi.mean <- rowGeoMeans(pig_monkey_mouse.rpms.cond1.a0.mouse_e55epi)

pig_monkey_mouse.rpms.cond1.a0.mean <- cbind(pig_monkey_mouse.rpms.cond1.a0.pig_e6icm.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.pig_e8epi.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.pig_e11epi.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.monkey_icm.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.monkey_pre_epi.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.monkey_poste_epi.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.mouse_e35icm.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.mouse_e45epi.mean,
                                             pig_monkey_mouse.rpms.cond1.a0.mouse_e55epi.mean)
colnames(pig_monkey_mouse.rpms.cond1.a0.mean) <- c("Pig_E6ICM","Pig_E8EPI","Pig_E11EPI","Monkey_ICM","Monkey_Pre-EPI","Monkey_PostE-EPI","Mouse_E3.5ICM","Mouse_E4.5EPI","Mouse_E5.5EPI")

rownames(pig_monkey_mouse.rpms.cond1.a0.mean) <- rownames(pigcell.geneid[pigcell.geneid[,1]%in%rownames(pig_monkey_mouse.rpms.cond1.a0.mean),])

library(corrplot)

pig_monkey_mouse.rpms.cond1.a0.corr <- cor(log2(pig_monkey_mouse.rpms.cond1.a0.mean+1))
corrplot(pig_monkey_mouse.rpms.cond1.a0.corr, method = "pie", hclust.method = "ward.D2", order="original", addrect = 3, cl.cex = 0.8, cl.lim = c(0,1), tl.cex = 1, tl.col = "black", mar = c(2,2,2,2), col = color.palette)

####
corr.input <- log2(pig_monkey_mouse.rpms.cond1.a0.mean[rownames(pig_monkey_mouse.rpms.cond1.a0.pc3sig),]+1)
corr.res <- cor(corr.input)
corr.test <- cor.mtest(corr.input)

corrplot(corr.res, method = "pie", hclust.method = "ward.D2", order="original", addrect = 3, cl.cex = 0.8, cl.lim = c(-1,1), tl.cex = 1, tl.col = "black", mar = c(2,2,2,2), col = color.palette)

####

layout(mat = matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE), widths = c(1,8), heights = c(2,8))

par(mai=c(0,0,0,0))
plot.new()

par(mar=c(0,0,0,0))
plot(c(1:11), rep(4,11), type = "n", frame = F, xaxt = "n", yaxt = "n", xlim = c(1,11), ylim = c(0,4))

rect(1, 0.5, 2, 1, col = "royalblue1", border = NA)
rect(2, 0.5, 3, 1, col = "lightgreen", border = NA)
rect(3, 0.5, 4, 1, col = "coral", border = NA)
rect(4, 0.5, 5, 1, col = "slateblue", border = NA)
rect(5, 0.5, 6, 1, col = "yellowgreen", border = NA)
rect(6, 0.5, 7, 1, col = "tomato", border = NA)
rect(7, 0.5, 8, 1, col = "skyblue", border = NA)
rect(8, 0.5, 9, 1, col = "seagreen1", border = NA)
rect(9, 0.5, 10, 1, col = "orange", border = NA)

text(1.5, 1.3, "ICM", font = 2, cex = 1.5)
text(2.5, 1.3, "EPI", font = 2, cex = 1.5)
text(3.5, 1.3, "EPI", font = 2, cex = 1.5)
text(4.5, 1.3, "ICM", font = 2, cex = 1.5)
text(5.5, 1.3, "EPI", font = 2, cex = 1.5)
text(6.5, 1.3, "EPI", font = 2, cex = 1.5)
text(7.5, 1.3, "ICM", font = 2, cex = 1.5)
text(8.5, 1.3, "EPI", font = 2, cex = 1.5)
text(9.5, 1.3, "EPI", font = 2, cex = 1.5)

text(1.5, 1.8, "", font = 2, cex = 1.5)
text(2.5, 1.8, "LB", font = 2, cex = 1.5)
text(3.5, 1.8, "Sph", font = 2, cex = 1.5)
text(4.5, 1.8, "", font = 2, cex = 1.5)
text(5.5, 1.8, "Pre-", font = 2, cex = 1.5)
text(6.5, 1.8, "PostE-", font = 2, cex = 1.5)
text(7.5, 1.8, "E3.5", font = 2, cex = 1.5)
text(8.5, 1.8, "E4.5", font = 2, cex = 1.5)
text(9.5, 1.8, "E5.5", font = 2, cex = 1.5)

arrows(1.1, 2.4, 3.9, 2.4, col = "navy", lwd = 4, length = 0.1, code = 3)
arrows(4.1, 2.4, 6.9, 2.4, col = "navy", lwd = 4, length = 0.1, code = 3)
arrows(7.1, 2.4, 9.9, 2.4, col = "navy", lwd = 4, length = 0.1, code = 3)

text(2.5, 3, "Pig", font = 2, cex = 1.5)
text(5.5, 3, "Monkey", font = 2, cex = 1.5)
text(8.5, 3, "Mouse", font = 2, cex = 1.5)

par(mai=c(0,0,0,0))
plot.new()

corrplot.mixed(corr.res,
               upper = "pie", upper.col = color.palette,
               lower = "circle", lower.col = color.palette,
               order="original",
               #p.mat = corr.test$p, sig.level = .05,
               cl.cex = 1.5, cl.lim = c(-1,1),
               tl.pos = "n", mar = c(0,0,0,0))

corrplot(corr.res, type = "lower", method = "number",
         col = "navy", add = T, diag = F,
         order="original", number.cex = 1.5,
         bg = NULL, addgrid.col = NA,
         cl.pos = "n", tl.pos = "n", mar = c(0,0,0,0))

