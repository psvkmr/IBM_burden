FYCO1 = "3: 45,917,899-45,995,824"
TBK1 = "12: 64,452,092-64,502,114"

fyco1.vars = filter(var.cleaned, `#CHROM` == 'chr3', (POS >= 45917899 & POS <= 45995824))
tbk1.vars = filter(var.cleaned, `#CHROM` == 'chr12', (POS >= 64452092 & POS <= 64502114))

fwrite(fyco1.vars, 'fyco1_variant_assoc_res.csv')
fwrite(tbk1.vars, 'tbk1_variant_assoc_res.csv')

tmp12.vep <- fread2('tmp12_var_filt.vcf.gz.vep.out.processed.csv')
tmp12.vep2 <- tmp12.vep[-1, ]

names(fyco1.vars)[1] <- 'CHROM'
names(tbk1.vars)[1] <- 'CHROM'

genes.df <- bind_rows(fyco1.vars, tbk1.vars)
info.df <- left_join(tmp12.vep, genes.df, by = c('CHROM', 'POS', 'REF', 'ALT'))
info.cleaned <- info.df[!is.na(info.df$ref_only_count),] %>% dplyr::select(-c(contains('LI'), count_hom, count_het, Phenotypes, panelWithQuality, ID.y, QUAL.y))

#fwrite(info.cleaned, 'fyco1_tbk1_fishers_with_anno.csv')

names(var.cleaned)[1] <- 'CHROM'

TBCD.vars = filter(var.cleaned, CHROM == 'chr17', (POS >= 82752065 & POS <= 82945914))
DLX6.vars = filter(var.cleaned, CHROM == 'chr7', (POS >= 97005553 & POS <= 97011040))
RNF225.vars = filter(var.cleaned, CHROM == 'chr19', (POS >= 58396090 & POS <= 58397079))
OR1S1.vars = filter(var.cleaned, CHROM == 'chr11', (POS >= 58212720 & POS <= 58216084))
SREBF1.vars = filter(var.cleaned, CHROM == 'chr17', (POS >= 17810399 & POS <= 17837011))
BCKDHB.vars = filter(var.cleaned, CHROM == 'chr6', (POS >= 80106647 & POS <= 80346270))
NECTIN2.vars = filter(var.cleaned, CHROM == 'chr19', (POS >= 44846175 & POS <= 44889223))
NRG2.vars = filter(var.cleaned, CHROM == 'chr5', (POS >= 139846779 & POS <= 140043299))
CARNMT1.vars = filter(var.cleaned, CHROM == 'chr9', (POS >= 74980790 & POS <= 75028423))
RUNX2.vars = filter(var.cleaned, CHROM == 'chr6', (POS >= 45328157 & POS <= 45664349))

gene.vars.list <- list(fyco1.vars, tbk1.vars, TBCD.vars, DLX6.vars, RNF225.vars, OR1S1.vars, SREBF1.vars, BCKDHB.vars, NECTIN2.vars, NRG2.vars, CARNMT1.vars, RUNX2.vars)

genes.df <- bind_rows(gene.vars.list)
genes.df$POS <- as.character(genes.df$POS)
info.df <- left_join(tmp12.vep, genes.df, by = c('CHROM', 'POS', 'REF', 'ALT'))
#info.cleaned <- info.df[!is.na(info.df$ref_only_count),] %>% dplyr::select(-c(count_hom, count_het, Phenotypes, panelWithQuality, ID.y, QUAL.y))
info.cleaned <- info.df[!is.na(info.df$ref_only_count),] %>% dplyr::select(-c(ID.y, QUAL.y))

#fwrite(info.cleaned, 'genes12_fishers_with_anno.csv')

# chr16	28499854	.	TCTC	T : LI1919 LI1986 LI2077 LI2101 LI2195 LI2199 LI2206  LI507  LI585  LI599
# chr16	28499876	.	TTCCTCC	TTCC : LI2087 LI2088  LI508  LI583  LI602
