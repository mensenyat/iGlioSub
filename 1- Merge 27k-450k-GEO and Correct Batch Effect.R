k27=read.table("GBM_Met_27k", sep=";", row.names=1, header=TRUE)
k450=read.table("GBM_Met_450k.csv", header=TRUE, row.names=1)
colnames(k450) = substr(colnames(k450), 1, 12)
k450 = k450[which(!colnames(k450)%in%colnames(k27))]
ClinAll=read.csv("Clinical_Data_TCGA.csv", sep=";", header=TRUE)

Clin27=subset(ClinAll, ClinAll$Patient.ID%in%colnames(k27))
Clin450=subset(ClinAll, ClinAll$Patient.ID%in%colnames(k450))

CLPer = nrow(Clin27[which(Clin27$Calculated.subtype=="Classical"),])/nrow(Clin27)
MEPer = nrow(Clin27[which(Clin27$Calculated.subtype=="Mesenchymal"),])/nrow(Clin27)
PNPer = nrow(Clin27[which(Clin27$Calculated.subtype=="Proneural"),])/nrow(Clin27)

c(CLPer, MEPer, PNPer)*nrow(Clin450)

CL450 = Clin450[which(Clin450$Calculated.subtype == "Classical"),]
CL450s = sample(CL450$Patient.ID, 10)

ME450 = Clin450[which(Clin450$Calculated.subtype == "Mesenchymal"),]
ME450s = sample(ME450$Patient.ID, 13)

PN450 = Clin450[which(Clin450$Calculated.subtype == "Proneural"),]
PN450s = sample(PN450$Patient.ID, 6)

Clin450sel = Clin450[which(Clin450$Patient.ID%in%c(as.character(CL450s), as.character(ME450s), as.character(PN450s))),]

k450sel = k450[which(colnames(k450)%in%Clin450sel$Patient.ID)]
k450sel = na.omit(k450sel)
  
GSE=read.table("GSE128654_All.csv", sep=";", header=TRUE, row.names=1)
GEOClin=read.csv("GEO122586.csv", row.names=1, sep=";", header=TRUE)
GEOClin = as.data.frame(as.matrix(t(GEOClin)))
table(GEOClin$Subtype)

CL450 = GEOClin[which(GEOClin$Subtype == "Classical"),]
CL450s = sample(rownames(CL450), 6)

ME450 = GEOClin[which(GEOClin$Subtype == "Mesenchymal"),]
ME450s = sample(rownames(ME450), 8)

PN450 = GEOClin[which(GEOClin$Subtype == "Proneural"),]
PN450s = sample(rownames(PN450), 4)

Clin450sel = GEOClin[which(rownames(GEOClin)%in%c(as.character(CL450s), as.character(ME450s), as.character(PN450s))),]

GSEsel = GSE[which(colnames(GSE)%in%rownames(Clin450sel))]
GSEsel = na.omit(GSEsel)

# Merge all datasets #
names450 = paste("k450", colnames(k450sel), sep="", collapse = NULL)
colnames(k450sel) = names450
names27 = paste("k27", colnames(k27), sep="", collapse = NULL)
colnames(k27) = names27
namesGEO = paste("GEO", colnames(GSEsel), sep="", collapse = NULL)
colnames(GSEsel) = namesGEO


k27450 = merge(k27, k450sel, by=0)
rownames(k27450) = k27450[,1]
k27450 = k27450[,-1]

data = merge(k27450, GSEsel, by=0)
rownames(data) = data[,1]
data = data[,-1]



library(sva)
cb.df.mdata <- cbind.data.frame("sample" = colnames(data), # exclude uid column, c(1)
                                "batch" = c(rep("k27", ncol(k27)), rep("k450", ncol(k450sel)), rep("GEO", ncol(GSEsel))))

cb.corr.model <- model.matrix(~1, data = cb.df.mdata)

rowna = rownames(data)
colna = colnames(data)
data=matrix(as.double(as.matrix(data)), ncol=ncol(data))
rownames(data) = rowna
colnames(data) = colna

cb.corr.counts = ComBat(dat=data,
                        batch=cb.df.mdata$batch,
                        mod=cb.corr.model,
                        par.prior=TRUE,
                        prior.plot=FALSE)
write.table(cb.corr.counts, file="Methylation_27_450_GEO.csv", sep=";")

# Plot tSNEs to identify changes due to batch #
tumorLabels=as.character(cb.df.mdata$batch)
plottingColors=c("Red", "cadetblue", "green")
names(plottingColors)=unique(tumorLabels)

library(Rtsne)
Tsne = Rtsne(as.matrix(t(cb.corr.counts)),perplexity=5)
# Static 2D plot of tSNE
pdf('Methylation_Corrected.pdf')
plot(Tsne$Y, main="Corrected Methylation", col=plottingColors[tumorLabels], pch=19)
legend('bottomright', legend=unique(tumorLabels), fill=plottingColors[unique(tumorLabels)], border=T, title='Subtype')
dev.off()

Tsne = Rtsne(as.matrix(t(data)),perplexity=5)
# Static 2D plot of tSNE
pdf('Methylation_Uncorrected.pdf')
plot(Tsne$Y, main="Uncorrected Methylation", col=plottingColors[tumorLabels], pch=19)
legend('bottomleft', legend=unique(tumorLabels), fill=plottingColors[unique(tumorLabels)], border=T, title='Subtype')
dev.off()


