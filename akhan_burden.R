library(tidyverse)
library(data.table)
library(readxl)

setwd("/array/psivakumar/other/akhan")

# remove contaminated samples

tmp1.samples <- fread("tmp1_samples.txt", header = F)
contamination.info <- fread("/array/psivakumar/191127_tracked.csv")
not.contaminated.samples <- filter(contamination.info, contamination < 0.03)$LOGID
not.contaminated.samples.in.tmp1 <- tmp1.samples[tmp1.samples$V1 %in% not.contaminated.samples, ]
#fwrite(not.contaminated.samples.in.tmp1, "notContaminatedToKeep.txt", col.names = F)

# new fam file

sample.ped <- fread("sample_ped.csv", header = F)[, 1:6]
id.to.log <- fread("unique_samples.csv")
id.to.log$V1 <- gsub("^", "LI", id.to.log$V1)
akhan.tmp1.fam <- fread("akhan_1_plink.fam")

new.fam <- left_join(akhan.tmp1.fam, id.to.log, by = c("V2" = "V1")) %>%
  left_join(sample.ped, by = c("V11" = "V2")) %>%
  dplyr::select("FID" = V11, "IID" = V2, "MATID" = V3, "FATID" = V4, "SEX" = V5, "PHENO" = V6) %>%
  arrange(IID)

new.fam$SEX <- gsub("other", "0", new.fam$SEX)

#fwrite(new.fam, "akhan_1_plink.fam", sep = '\t', col.names = F)

### ISTATS

istats <- fread("pre_QC_istats.txt")

detectOutliers <- function(x) {
  mean.x <- mean(x)
  sd.x <- sd(x)
  outliers <- c()
  for (i in 1:length(x)){
    if (x[i] > mean.x + 3*sd.x) {
      outliers <- c(outliers, x[i])
    } else if (x[i] < mean.x - 3*sd.x) {
      outliers <- c(outliers, x[i])
    }
  }
  outliers.rows <- match(outliers, x)
  return(outliers.rows)
}

filter.columns <- c(2:7)
istats.outliers <- lapply(istats[, ..filter.columns], detectOutliers)
outlier.samples <- lapply(istats.outliers, function(x) istats[unlist(x), 1])
all.outlier.samples <- unique(unlist(outlier.samples))
outlier.samples.df <- filter(new.fam, IID %in% all.outlier.samples)[, 1:2]
#fwrite(outlier.samples.df, "istats_outliers_to_remove.txt", col.names = FALSE, sep = '\t')

# high missing samples
ind.miss <- fread("individual_missing_rate.imiss")

## het outliers

hets <- fread("het_statistics.het", header = T)
hets$HET_RATE = (hets$`N(NM)` - hets$`O(HOM)`) / hets$`N(NM)`
het.fail = subset(hets, (hets$HET_RATE < mean(hets$HET_RATE)-3*sd(hets$HET_RATE)) | (hets$HET_RATE > mean(hets$HET_RATE)+3*sd(hets$HET_RATE)))
het.fail$HET_DST = (het.fail$HET_RATE - mean(hets$HET_RATE)) / sd(hets$HET_RATE)
het.fail <- het.fail[, 1:2]
#write.table(het.fail, "hets_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


## relatives

rel.check <- fread("relatives_check.genome")
duplicate.samples <- filter(rel.check, PI_HAT > 0.95)$IID1
rel.to.remove <- filter(rel.check, PI_HAT > 0.25)[, 1:2] %>% unique.data.frame()
#write.table(rel.to.remove, "related_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


## sex check

sex.check <- read.table("sex_check.sexcheck", header = T)
# imputed sexes all look fine so taken no further


## pop stratification

pop.strat <- read.csv("peddy_emer.het_check.csv")
pop.dist <- table(pop.strat$ancestry.prediction)
european <- filter(pop.strat, ancestry.prediction == "EUR")$sample_id
#write.table(european, "european_to_keep.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


### PCA

eigenvecs <- read.table("pca.eigenvec")
eigenvals <- read.table("pca.eigenval")

sum.eig <- sum(eigenvals$V1)

sum.eigs <- lapply(eigenvals$V1, function(x){
  rt<-(x/sum.eig)*100
  rt<-round(rt)
  return(rt)
})

pca.plot <- ggplot(eigenvecs, aes(V3, V4)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",sum.eigs[[1]],"% variance")) +
  ylab(paste0("PC2: ",sum.eigs[[2]],"% variance"))

pca.filter.columns <- c(3:12)
pc.outliers <- lapply(eigenvecs[, pca.filter.columns], detectOutliers)
pca.outlier.samples <- lapply(pc.outliers, function(x) eigenvecs[unlist(x), 2])
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
#pca.outlier.samples.df <- data.frame(fid = 0, iid = pca.all.outlier.samples)
pca.outlier.samples.df <- data.frame()
#write.table(pca.outlier.samples.df, "pca_samples_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


# create separate case and control samples list

all.samples.final <- read.table("samples_final_set.txt")
case.samples <- filter(all.samples.final, V1 %in% filter(new.fam, PHENO == 2)$IID)
control.samples <- filter(all.samples.final, V1 %in% filter(new.fam, PHENO == 1)$IID)
#write.table(case.samples, "case_samples_final_set.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(control.samples, "control_samples_final_set.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


## pre vs post QC PCA

#pre_qc_pca <- left_join(pop.strat, new.fam, by = c("sample_id" = "V2")) %>%
#  ggplot(aes(PC1, PC2, colour = as.factor(V6))) + geom_point()
#post_qc_pca <- left_join(pop.strat, new.fam, by = c("sample_id" = "V2")) %>%
#  filter(sample_id %in% all.samples.final$V1) %>%
#  ggplot(aes(PC1, PC2, colour = as.factor(V6))) + geom_point()

### create list of phenotype specific low quality variants to remove

cases.depth <- read.table("case_new_depth.ldepth", header = T)
cases.depth.pass <- filter(cases.depth, (SUM_DEPTH / nrow(case.samples)) > 10)
controls.depth <- read.table("control_new_depth.ldepth", header = T)
controls.depth.pass <- filter(controls.depth, (SUM_DEPTH / nrow(control.samples)) > 10)
depth.pass <- bind_rows(cases.depth.pass[, 1:2], controls.depth.pass[, 1:2]) %>% unique.data.frame()
#write.table(depth.pass, "read_depth_pass_to_include.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

### create refflat and set files for rvtests burden

library(biomaRt)

gene.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
gene.dataset <- useDataset(mart = gene.mart, dataset = "hsapiens_gene_ensembl")
genes.pos <- getBM(mart = gene.dataset,
                   attributes = c("external_gene_name","chromosome_name", "start_position", "end_position", "gene_biotype"),
                   filters = c("biotype"),
                   values = list("protein_coding"))
tmp.set.file <- filter(genes.pos, chromosome_name %in% seq(1, 22, 1))
tmp.set.file$chromosome_name <- gsub("^", "chr", tmp.set.file$chromosome_name)
set.file <- unite(tmp.set.file, "tmp", c("chromosome_name", "start_position"), sep = ":") %>%
  unite("pos", c("tmp", "end_position"), sep = "-")

#write.table(set.file, "hg38_genes.set", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

chrs <- seq(1, 22, 1) %>% paste('chr', ., sep = '')
#refflat <- read.table("J:/NGS_Reference/refFLAT/refFlat_hg38.txt")
#refflat22 <- filter(refflat, V3 %in% chrs)

#write.table(refflat22, "hg38_refFlat_chr1to22.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


library(bigreadr)

nonsyn.vcf <- fread2("rare_nonsyn.vcf", sep = '\t', skip = 3444)
syn.vcf <- fread2("rare_syn.vcf", sep = '\t', skip = 3443)
full.vcf <- fread2('tmp12_var_filt.vcf', sep = '\t', skip = 3425)

codeGenotypes <- function(x) {
  y <- strsplit(x, split = ":")[[1]][1]
  if (y == "./.") {
    return(0)
  } else if (y == "0/0") {
    return(0)
  } else if (y == "0/1") {
    return(1)
  } else if (y == "1/0") {
    return(1)
  } else if (y == "1/1") {
    return(2)
  } else {
    return(NA)
  }
}

#test.mat <- apply(test.vcf[, 10:ncol(test.vcf)], c(1, 2), codeGenotypes)
nonsyn.mat <- apply(nonsyn.vcf[, 10:ncol(nonsyn.vcf)], c(1, 2), codeGenotypes)
syn.mat <- apply(syn.vcf[, 10:ncol(syn.vcf)], c(1, 2), codeGenotypes)
full.mat <- apply(full.vcf[, 10:ncol(full.vcf)], c(1, 2), codeGenotypes)

varAddGenes <- function(vcf, setfile) {
  gene.mat.list <- list(length = nrow(setfile))
  for (i in 1:nrow(setfile)) {
    gene.mat <- which(vcf$`#CHROM` == setfile[i, 2] & vcf$POS > setfile[i, 3] & vcf$POS < setfile[i, 4])
    gene.mat.list[[i]] <- gene.mat
  }
  return(gene.mat.list)
}

#test.gene.list <- varAddGenes(test.vcf, tmp.set.file)
nonsyn.gene.list <- varAddGenes(nonsyn.vcf, tmp.set.file)
syn.gene.list <- varAddGenes(syn.vcf, tmp.set.file)
full.gene.list <- varAddGenes(full.vcf, tmp.set.file)

#getAltCounts <- function(mat, genelist) {
#  gene.data <- data.frame()
#  case.mat <- mat[, unlist(case.samples)]
#  control.mat <- mat[, unlist(control.samples)]
#  for (i in 1:length(genelist)) {
#    nvar <- length(genelist[[i]])
#    alt.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] != 0))
#    het.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] == 1))
#    hom.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] == 2))
#    control.alt <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] != 0))
#    control.het <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] == 1))
#    control.hom <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] == 2))
#    case.alt <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] != 0))
#    case.het <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] == 1))
#    case.hom <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] == 2))
#    gene.df <- c(nvar, alt.count, het.count, hom.count, control.alt, control.het, control.hom, case.alt, case.het, case.hom)
#    gene.data <- rbind(gene.data, gene.df)
#  }
#  names(gene.data) <- c("nvar.count", "alt.count", "het.count", "hom.count", "control.alt.count", "control.het.count", "control.hom.count", "case.alt.count", "case.het.count", "case.hom.count")
#  #mutate(gene.data, ref.only.count = nrow(all.samples.final) - alt.count, control.ref.only.count = nrow(control.samples) - control.alt.count, case.ref.only.count = nrow(case.samples) - case.alt.count)
#  #dplyr::select(gene.data, c(1, 11, 2:4, 12, 5:7, 13, 8:10))
#  return(gene.data)
#}

getAltCounts <- function(mat, genelist) {
  gene.data <- data.frame()
  case.mat <- mat[, unlist(case.samples)]
  control.mat <- mat[, unlist(control.samples)]
  for (i in 1:length(genelist)) {
    nvar <- length(genelist[[i]])
    ref.only.count <- sum(colSums(mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    has.alt.count <- sum(colSums(mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    control.ref.only.count <- sum(colSums(control.mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    control.has.alt.count <- sum(colSums(control.mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    case.ref.only.count <- sum(colSums(case.mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    case.has.alt.count <- sum(colSums(case.mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    het.count <- length(which(mat[genelist[[i]], , drop = FALSE] == 1))
    hom.count <- length(which(mat[genelist[[i]], , drop = FALSE] == 2))
    control.het.count <- length(which(control.mat[genelist[[i]], , drop = FALSE] == 1))
    control.hom.count <- length(which(control.mat[genelist[[i]], , drop = FALSE] == 2))
    case.het.count <- length(which(case.mat[genelist[[i]], , drop = FALSE] == 1))
    case.hom.count <- length(which(case.mat[genelist[[i]], , drop = FALSE] == 2))
    gene.df <- c(nvar, ref.only.count, has.alt.count, control.ref.only.count, control.has.alt.count, case.ref.only.count, case.has.alt.count, het.count, hom.count, control.het.count, control.hom.count, case.het.count, case.hom.count)
    gene.data <- rbind(gene.data, gene.df)
  }
  names(gene.data) <- c("nvar_count", "ref_only_count", "alt_count", "control_ref_only_count", "control_alt_count", "case_ref_only_count", "case_alt_count", "het_count", "hom_count", "control_het_count", "control_hom_count", "case_het_count", "case_hom_count")
  return(gene.data)
}

#test.df <- getAltCounts(test.mat, test.gene.list)
nonsyn.df <- getAltCounts(nonsyn.mat, nonsyn.gene.list)
syn.df <- getAltCounts(syn.mat, syn.gene.list)
full.df <- getAltCounts(full.mat, full.gene.list)

#test.df.join <- cbind(tmp.set.file, test.df)
nonsyn.df.join <- cbind(tmp.set.file, nonsyn.df)
syn.df.join <- cbind(tmp.set.file, syn.df)
full.df.join <- cbind(tmp.set.file, full.df)

#test.assoc.join <- left_join(test.df.join, nonsyn.assoc, by = c("external_gene_name" = "Range"))
#nonsyn.assoc.join <- left_join(nonsyn.df.join, nonsyn.assoc, by = c("external_gene_name" = "Range"))
#syn.assoc.join <- left_join(syn.df.join, syn.assoc, by = c("external_gene_name" = "Range"))

runFishers <- function(df_join) {
  fishers.res <- data.frame()
  for (i in 1:nrow(df_join)) {
    test.res <- fisher.test(matrix(c(df_join[i, 9], df_join[i, 10], df_join[i, 11], df_join[i, 12]), nrow = 2, byrow = T))
    res.vec <- c(test.res[[1]], test.res[[3]], test.res[[2]][1], test.res[[2]][2])
    fishers.res <- rbind(fishers.res, res.vec)
  }
  names(fishers.res) <- c("pvalue", "or", "lower_ci", "upper_ci")
  fishers.res <- cbind(df_join, fishers.res)
  return(fishers.res)
}

#test.fishers.res <- runFishers(test.df.join)
nonsyn.fishers.res <- runFishers(nonsyn.df.join)
syn.fishers.res <- runFishers(syn.df.join)
full.fishers.res <- runFishers(full.df.join)

# genes with duplicate location regions
#problem.genes <- filter(nonsyn.assoc.join, paste0(chromosome_name, ":", start_position, "-", end_position) != RANGE)$external_gene_name %>% unique()

nonsyn.cleaned <- nonsyn.fishers.res %>% #[!nonsyn.fishers.res$external_gene_name %in% problem.genes, ] %>%
  filter(alt_count != 0) %>%
  #select(c(1:11, 18:19, 23:26)) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
  arrange(pvalue)

#fwrite(nonsyn.cleaned, "rare_nonsynonymous_burden_results.csv")

syn.cleaned <- syn.fishers.res %>% #[!syn.fishers.res$external_gene_name %in% problem.genes, ] %>%
  filter(alt_count != 0) %>%
  #select(c(1:11, 18:19, 23:26)) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
  arrange(pvalue)

#fwrite(syn.cleaned, "rare_synonymous_burden_results.csv")

getVarAltCounts <- function(mat) {
  var.data <- data.frame()
  case.mat <- mat[, case.samples$V1]
  control.mat <- mat[, control.samples$V1]
  for (i in 1:dim(mat)[1]) {
    ref.only.count <- sum(colSums(mat[i, , drop = FALSE] != 0) == 0)
    has.alt.count <- sum(colSums(mat[i, , drop = FALSE] != 0) != 0)
    control.ref.only.count <- sum(colSums(control.mat[i, , drop = FALSE] != 0) == 0)
    control.has.alt.count <- sum(colSums(control.mat[i, , drop = FALSE] != 0) != 0)
    case.ref.only.count <- sum(colSums(case.mat[i, , drop = FALSE] != 0) == 0)
    case.has.alt.count <- sum(colSums(case.mat[i, , drop = FALSE] != 0) != 0)
    het.count <- length(which(mat[i, , drop = FALSE] == 1))
    hom.count <- length(which(mat[i, , drop = FALSE] == 2))
    control.het.count <- length(which(control.mat[i, , drop = FALSE] == 1))
    control.hom.count <- length(which(control.mat[i, , drop = FALSE] == 2))
    case.het.count <- length(which(case.mat[i, , drop = FALSE] == 1))
    case.hom.count <- length(which(case.mat[i, , drop = FALSE] == 2))
    var.df <- c(ref.only.count, has.alt.count, control.ref.only.count, control.has.alt.count, case.ref.only.count, case.has.alt.count, het.count, hom.count, control.het.count, control.hom.count, case.het.count, case.hom.count)
    var.data <- rbind(var.data, var.df)
  }
  names(var.data) <- c("ref_only_count", "alt_count", "control_ref_only_count", "control_alt_count", "case_ref_only_count", "case_alt_count", "het_count", "hom_count", "control_het_count", "control_hom_count", "case_het_count", "case_hom_count")
  return(var.data)
}

var.mat <- apply(full.vcf[, 10:ncol(full.vcf)], c(1, 2), codeGenotypes)
var.df <- getVarAltCounts(var.mat)
var.df.join <- cbind(full.vcf[, 1:6], var.df)
var.fishers.res <- runFishers(var.df.join)
var.cleaned <- var.fishers.res %>%
  filter(alt_count != 0) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
  arrange(pvalue)
fwrite(var.cleaned, "single_variant_assoc_results.csv")


## qqplots

library(ggrepel)

qqplotter <- function(df) {
  df <- df[!is.na(df$pvalue), ]
  df <- mutate(df, exp.pvalues = (rank(df$pvalue, ties.method="first")+.5)/(length(df$pvalue)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(pvalue), label = external_gene_name)) +
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") +
    geom_abline(slope = 1, intercept = 0) +
    geom_label_repel(data = arrange(df, pvalue)[1:5, ], box.padding = 0.5) +
    theme_classic()
}


nonsyn.qq <- qqplotter(nonsyn.cleaned)
syn.qq <- qqplotter(syn.cleaned)

nonsyn.chisq <- qchisq(1-nonsyn.cleaned$pvalue, 1)
nonsyn.lambda <- median(nonsyn.chisq) / qchisq(0.5,1)

syn.chisq <- qchisq(1-syn.cleaned$pvalue, 1)
syn.lambda <- median(syn.chisq) / qchisq(0.5,1)

# go terms

library(gProfileR)

goTerms <- function(df_cleaned) {
  go.terms <- gprofiler(query = df_cleaned[1:2000, 1],
                        organism = "hsapiens",
                        underrep = F,
                        min_set_size = 5,
                        custom_bg = df_cleaned$external_gene_name)
}

nonsyn.go <- goTerms(nonsyn.cleaned)
syn.go <- goTerms(syn.cleaned)
