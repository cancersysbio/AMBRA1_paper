# figure 5 panel h
# The AMBRA1 E3 ligase adaptor regulates Cyclin D protein stability
# Chaikovsky et al. Nature 2021
# contact: preethi.srinivasan@stanford.edu 

library(MASS)
library(gplots)
library(RColorBrewer)




load("AMBRA1_expression_analysis.RData")



clinic.kras = clinic.tmp[which(clinic.tmp$KRAS_G12_status==1),]
clinic.egfr = clinic.tmp[which(clinic.tmp$EGFR_hotspot_status==1),]
clinic.kraswt = clinic.tmp[which(clinic.tmp$KRAS_status==0),]

dim(clinic.kras)
dim(clinic.egfr)
dim(clinic.kraswt)


#to generate models for KRAS wildtype and EGFR mutant, substitute clinic.kras with clinic.kraswt and clinic.egfr below respectively.
z= intersect(substr(colnames(exp_luad),1,12),clinic.kras$bcr_patient_barcode)
clinic.kras = clinic.kras[match(z,clinic.kras$bcr_patient_barcode),]
dim(clinic.kras )
rppa_luad = rppa_luad[,match(z,substr(colnames(rppa_luad),1,12))]
exp_luad = exp_luad[,match(z,substr(colnames(exp_luad),1,12))]
samplesNUM_1 = samplesNUM_1[match(z,substr(rownames(samplesNUM_1),1,12)),,drop=F]

dim(exp_luad)
dim(rppa_luad)
dim(samplesNUM_1)

cn_adjust <- lm(as.numeric(rppa_luad["Cyclin_D1-R-V",]) ~ as.numeric(samplesNUM_1[,"595"]),)
cn.resid <- residuals(cn_adjust)


df_exp <-data.frame(exp_luad)
df_exp <- as.data.frame(t(df_exp))
pathway_genes = c("ENSG00000110497", "ENSG00000147889", "ENSG00000147883", "ENSG00000123080", 
                  "ENSG00000118971", "ENSG00000112576", "ENSG00000136997", "ENSG00000135446", 
                  "ENSG00000105810", "ENSG00000110092", "ENSG00000139687")

df_exp<-df_exp[complete.cases(df_exp),]
df_exp<-df_exp[, pathway_genes]

#get model for KRAS
full.model <- lm(cn.resid ~., data = df_exp)
step.model <- stepAIC(full.model, direction = "both", 
                      trace = FALSE)
summary(step.model)


fit <- lm(cn.resid ~ 
            as.numeric(exp_luad["ENSG00000110497",]) +
            as.numeric(exp_luad["ENSG00000112576",]) +
            as.numeric(exp_luad["ENSG00000136997",]) +
            as.numeric(exp_luad["ENSG00000110092",]))

summary(fit)

##Plot heatmap KRAS mutant
df_exp_plot<-df_exp[, c("ENSG00000110497", "ENSG00000112576", "ENSG00000136997", "ENSG00000110092")]
protein_cyclin <- rppa_luad["Cyclin_D1-R-V",]
protein_cyclin <- as.data.frame(t(protein_cyclin))
rownames(protein_cyclin) <-substr(rownames(protein_cyclin),1,12)
rownames(df_exp_plot) <-substr(rownames(df_exp_plot),1,12)
rownames(df_exp_plot) <- gsub('\\.', '-', rownames(df_exp_plot))
df_exp_plot<-merge(df_exp_plot, protein_cyclin, by="row.names", all=TRUE)
df_exp_plot<-subset(df_exp_plot, select=c( "ENSG00000110497", "ENSG00000112576", "ENSG00000136997", "ENSG00000110092", "Cyclin_D1-R-V"))

df_exp_plot <- df_exp_plot[order(df_exp_plot$`Cyclin_D1-R-V`),]#, -df_exp_plot$ENSG00000110497),]
colnames(df_exp_plot)<-c('AMBRA1', 'CCND3', 'MYC', 'CCND1', 'Cyclin_D1')
df_exp_plot<-subset(df_exp_plot, select=c("Cyclin_D1", 'AMBRA1', 'CCND1', 'MYC', 'CCND3'))


heatmap.2(t(as.matrix(df_exp_plot)), Colv = NA, Rowv = NA, col=brewer.pal(9,"Blues"), scale="row",symbreaks = FALSE,labCol=FALSE, 
          density.info="none", trace="none")


##Model and plot heatmap EGFR mutant
fit <- lm(cn.resid ~ #EGFR mut
            as.numeric(exp_luad["ENSG00000147889",]) +
            as.numeric(exp_luad["ENSG00000147883",]) +
            as.numeric(exp_luad["ENSG00000123080",]) +
            as.numeric(exp_luad["ENSG00000135446",]) +
            as.numeric(exp_luad["ENSG00000110092",]) +
            as.numeric(exp_luad["ENSG00000110497",]))

summary(fit)



df_exp_plot<-df_exp[, c("ENSG00000110497", "ENSG00000110092", "ENSG00000135446", "ENSG00000123080", "ENSG00000147883", "ENSG00000147889")]
protein_cyclin <- rppa_luad["Cyclin_D1-R-V",]
protein_cyclin <- as.data.frame(t(protein_cyclin))
rownames(protein_cyclin) <-substr(rownames(protein_cyclin),1,12)
rownames(df_exp_plot) <-substr(rownames(df_exp_plot),1,12)
rownames(df_exp_plot) <- gsub('\\.', '-', rownames(df_exp_plot))
df_exp_plot<-merge(df_exp_plot, protein_cyclin, by="row.names", all=TRUE)
df_exp_plot<-subset(df_exp_plot, select=c("ENSG00000110497", "ENSG00000110092", "ENSG00000135446", "ENSG00000123080", "ENSG00000147883", "ENSG00000147889", "Cyclin_D1-R-V"))

df_exp_plot <- df_exp_plot[order(df_exp_plot$`Cyclin_D1-R-V`),]
colnames(df_exp_plot)<-c('AMBRA1', 'CCND1', 'CDK4', 'CDKN2C', 'CDKN2B', 'CDKN2A', 'Cyclin_D1')
df_exp_plot<-subset(df_exp_plot, select=c("Cyclin_D1", 'AMBRA1', 'CCND1', 'CDK4', 'CDKN2C', 'CDKN2B', 'CDKN2A'))


heatmap.2(t(as.matrix(df_exp_plot)), Colv = NA, Rowv = NA, col=brewer.pal(9,"Blues"), scale="row",symbreaks = FALSE,labCol=FALSE, 
          density.info="none", trace="none")


