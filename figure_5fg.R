# figure 5 panel f, g
# The AMBRA1 E3 ligase adaptor regulates Cyclin D protein stability
# Chaikovsky et al. Nature 2021
# contact: seoane@stanford.edu

# LUAD TCGA data
# clinical data obtained from 29625055 (Liu et al. Cell 2018)
# AMBRA1_z log2 FPKM-UQ rnaseq data from gdc.cancer.gov
# AMBRA1_CN AMBRA1 data from segmented file from gdc.cancer.gov
# AMBRA1_adj copy number adjusted expression data (lm(AMBRA1_z~AMBRA1_CN,clinic)$residuals)
# AMBRA1_res1_B 1/3 high and 1/3 low quantiles from AMBRA1_adj
# AMBRA1_res2_B 1/3 high and 1/3 low quantiles from AMBRA1_z
# KRAS_mut KRAS MC3 mutation calls from http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
# KRAS_mut_ms KRAS missense
# KRAS_mut_ms_G12 KRAS missesne mutations in p.G12A, p.G12C, p.G12D, p.G12S and p.G12v

library(survival)
library(survminer)

clinic = readRDS("clinic_molecular_luad_dataframe.Rds")

low_exp = "#66D8EE"
med_exp = "#8FB1DD"
high_exp= "#488CCA"

clinic.kras = clinic[which(clinic$KRAS_mut_ms_G12),]
qt = quantile(clinic.kras$AMBRA1_adj,probs = c(1/3,2/3),na.rm = T)
clinic.kras$AMBRA1_res1_B = NA
clinic.kras$AMBRA1_res1_B[clinic.kras$AMBRA1_adj>qt[2]]="High"
clinic.kras$AMBRA1_res1_B[clinic.kras$AMBRA1_adj<qt[1]]="Low"

clinic.kras$AGE=clinic.kras$age_at_initial_pathologic_diagnosis
clinic.kras$STAGE = factor(clinic.kras$PATH_STAGE)
clinic.kras$GENDER = factor(clinic.kras$gender)

m0.krasmut<- survfit(Surv(os.t.10y, os.e.10y) ~ AMBRA1_res1_B, data = clinic.kras)

m1.kras.coxph = coxph(Surv(os.t.10y,os.e.10y==1)~AMBRA1_adj+AGE+STAGE+GENDER,clinic.kras,x=T,y=T)

#Figure 5f
ggsurvplot(m0.krasmut,pval = T,risk.table = T,palette = c(high_exp,low_exp),risk.table.caption="KRAS mut",
           legend.title="AMBRA1",legend.labs=c("high","low"),
           font.caption=12,font.x=10,font.y=10,fonts.tickslab=10,pval.size=5,risk.table.fontsize=4,tables.height=0.3)

#figure 5g
ggforest(m1.kras.coxph,fontsize = 0.7)
