library(tidyverse)

# download Supplementary Table 4 from BB Sun, et al. Genomic atlas of the human plasma proteome. Nature 2018 ( https://www.nature.com/articles/s41586-018-0175-2#Sec37 ) and store as tab-delimited text file “st4.txt”.

st4_1 <- read_delim(file = "path_to_directory/st4.txt", delim = "\t", col_names = T)

st4_2 <- st4_1 %>% add_column(z = abs(st4_1$beta/st4_1$SE))

st4_3 <- st4_2 %>% arrange(desc(z))

st4_4 <- st4_3 %>% distinct(UniProt, .keep_all = TRUE)

# download and unzip the file “Phelan_Archive.zip” from ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/PhelanCM_28346442_GCST004462 .  I use a version of these data where the data for each chromosome is stored as a separate file named in the format Summary_chr*.txt where * is the chromosome number from 1 to 22 (chromosome 23 is not used in this analysis since the file “st4.txt” above only contains data for chromosomes 1 to 22).  Replace the number 1 with the numbers 2 to 22 and iterate the next two lines of code to obtain tibbles chr1x to chr22x.

chr1 <- read_delim(file = " path_to_directory/Summary_chr1.txt", delim = ",", col_names = T)

chr1x <- inner_join(st4_4, chr1, by = c("Chromosome","Position"))

d <- bind_rows(chr1x,chr2x,chr3x,chr4x,chr5x,chr6x,chr7x,chr8x,chr9x,chr10x,chr11x,chr12x,chr13x,chr14x,chr15x,chr16x,chr17x,chr18x,chr19x,chr20x,chr21x,chr22x)

# for tibble d, please see file interval_ocac.txt in this repository.

d1 <- d %>% filter(EAF.y > 0.01)

d2 <- d1 %>% filter(R2_oncoarray > 0.8)

d3 <- d2 %>% filter(INFO > 0.8)

d4 <- d3 %>% filter(EAF.x > 0.01)

d5 <- d4 %>% filter(abs(EAF.x-EAF.y)<0.05)

d6 <- d5 %>% filter(`Effect Allele (EA)`==Effect)

d6 <- d6 %>% filter(`Other Allele (OA)`==Baseline)

d6 <- add_column(d6, all_beta=d6$overall_OR/d6$beta)

d6 <- add_column(d6, all_se=d6$overall_SE/abs(d6$beta))

d6 <- add_column(d6, all_p=pnorm(abs(d6$all_beta)/d6$all_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, all_fdr=p.adjust(d6$all_p,method="BH"))

d6 <- add_column(d6, all_or=exp(d6$all_beta))

d6 <- add_column(d6, all_lcl=exp(d6$all_beta-(1.96*d6$all_se)))

d6 <- add_column(d6, all_ucl=exp(d6$all_beta+(1.96*d6$all_se)))

d6 <- add_column(d6, hgsoc_beta=d6$serous_hg_OR/d6$beta)

d6 <- add_column(d6, hgsoc_se=d6$serous_hg_SE/abs(d6$beta))

d6 <- add_column(d6, hgsoc_p=pnorm(abs(d6$hgsoc_beta)/d6$hgsoc_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, hgsoc_fdr=p.adjust(d6$hgsoc_p,method="BH"))

d6 <- add_column(d6, hgsoc_or=exp(d6$hgsoc_beta))

d6 <- add_column(d6, hgsoc_lcl=exp(d6$hgsoc_beta-(1.96*d6$hgsoc_se)))

d6 <- add_column(d6, hgsoc_ucl=exp(d6$hgsoc_beta+(1.96*d6$hgsoc_se)))

d6 <- add_column(d6, lgsoc_beta=d6$serouslowgrade_OR/d6$beta)

d6 <- add_column(d6, lgsoc_se=d6$serouslowgrade_SE/abs(d6$beta))

d6 <- add_column(d6, lgsoc_p=pnorm(abs(d6$lgsoc_beta)/d6$lgsoc_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, lgsoc_fdr=p.adjust(d6$lgsoc_p,method="BH"))

d6 <- add_column(d6, lgsoc_or=exp(d6$lgsoc_beta))

d6 <- add_column(d6, lgsoc_lcl=exp(d6$lgsoc_beta-(1.96*d6$lgsoc_se)))

d6 <- add_column(d6, lgsoc_ucl=exp(d6$lgsoc_beta+(1.96*d6$lgsoc_se)))

d6 <- add_column(d6, lmpsoc_beta=d6$serous_lmp_OR/d6$beta)

d6 <- add_column(d6, lmpsoc_se=d6$serous_lmp_SE/abs(d6$beta))

d6 <- add_column(d6, lmpsoc_p=pnorm(abs(d6$lmpsoc_beta)/d6$lmpsoc_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, lmpsoc_fdr=p.adjust(d6$lmpsoc_p,method="BH"))

d6 <- add_column(d6, lmpsoc_or=exp(d6$lmpsoc_beta))

d6 <- add_column(d6, lmpsoc_lcl=exp(d6$lmpsoc_beta-(1.96*d6$lmpsoc_se)))

d6 <- add_column(d6, lmpsoc_ucl=exp(d6$lmpsoc_beta+(1.96*d6$lmpsoc_se)))

d6 <- add_column(d6, mucinous_beta=d6$mucinous_OR/d6$beta)

d6 <- add_column(d6, mucinous_se=d6$mucinous_SE/abs(d6$beta))

d6 <- add_column(d6, mucinous_p=pnorm(abs(d6$mucinous_beta)/d6$mucinous_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, mucinous_fdr=p.adjust(d6$mucinous_p,method="BH"))

d6 <- add_column(d6, mucinous_or=exp(d6$mucinous_beta))

d6 <- add_column(d6, mucinous_lcl=exp(d6$mucinous_beta-(1.96*d6$mucinous_se)))

d6 <- add_column(d6, mucinous_ucl=exp(d6$mucinous_beta+(1.96*d6$mucinous_se)))

d6 <- add_column(d6, lmpmuc_beta=d6$mucinous_lmp_OR/d6$beta)

d6 <- add_column(d6, lmpmuc_se=d6$mucinous_lmp_SE/abs(d6$beta))

d6 <- add_column(d6, lmpmuc_p=pnorm(abs(d6$lmpmuc_beta)/d6$lmpmuc_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, lmpmuc_fdr=p.adjust(d6$lmpmuc_p,method="BH"))

d6 <- add_column(d6, lmpmuc_or=exp(d6$lmpmuc_beta))

d6 <- add_column(d6, lmpmuc_lcl=exp(d6$lmpmuc_beta-(1.96*d6$lmpmuc_se)))

d6 <- add_column(d6, lmpmuc_ucl=exp(d6$lmpmuc_beta+(1.96*d6$lmpmuc_se)))

d6 <- add_column(d6, cc_beta=d6$clearcell_OR/d6$beta)

d6 <- add_column(d6, cc_se=d6$clearcell_SE/abs(d6$beta))

d6 <- add_column(d6, cc_p=pnorm(abs(d6$cc_beta)/d6$cc_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, cc_fdr=p.adjust(d6$cc_p,method="BH"))

d6 <- add_column(d6, cc_or=exp(d6$cc_beta))

d6 <- add_column(d6, cc_lcl=exp(d6$cc_beta-(1.96*d6$cc_se)))

d6 <- add_column(d6, cc_ucl=exp(d6$cc_beta+(1.96*d6$cc_se)))

d6 <- add_column(d6, endo_beta=d6$endometrioid_OR/d6$beta)

d6 <- add_column(d6, endo_se=d6$endometrioid_SE/abs(d6$beta))

d6 <- add_column(d6, endo_p=pnorm(abs(d6$endo_beta)/d6$endo_se, lower.tail=FALSE) * 2)

d6 <- add_column(d6, endo_fdr=p.adjust(d6$endo_p,method="BH"))

d6 <- add_column(d6, endo_or=exp(d6$endo_beta))

d6 <- add_column(d6, endo_lcl=exp(d6$endo_beta-(1.96*d6$endo_se)))

d6 <- add_column(d6, endo_ucl=exp(d6$endo_beta+(1.96*d6$endo_se)))
