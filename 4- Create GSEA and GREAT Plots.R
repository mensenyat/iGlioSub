library(ggplot2)

setwd("D:/R/GBM-LGG/GBM Subtypes/GSEA")
CL = read.table("Classical - Pathways to plot.csv", sep=";", header = TRUE) # Table obtained from metascape containing only the pathways that we want to plot
CL$N = c(1:nrow(CL))

p = ggplot(CL, aes(x=Enrichment, y=N, color = LogP)) + # Enrichment is the LogP of Up or Down. In case of being downreglated is converted to negative
  geom_point(size=(CL$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(-max(abs(CL$Enrichment)), max(abs(CL$Enrichment))) +
  theme_minimal(base_size = 20)
p
ggsave(plot=p,height=NA,width=NA, filename="GSEA_Classical.pdf", useDingbats=FALSE)

# Plot to create manually a legend for the Figure
df = data.frame(Enrichment=1, ev2=c(40,80,120), ev3=1:3)
pd = ggplot(df, aes(x=Enrichment, y=ev3, color = -10)) + 
  geom_point(size=(df$ev2*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  theme_minimal()
ggsave(plot=pd,height=NA,width=NA, filename="Points_for_legend.pdf", useDingbats=FALSE) 


#### MESENCHYMAL ####
ME = read.table("Mesenchymal - Pathways to plot.csv", sep=";", header = TRUE)
ME$N = c(1:nrow(ME))

p = ggplot(ME, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(ME$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(-max(abs(ME$Enrichment)), max(abs(ME$Enrichment))) +
  geom_vline(xintercept = 0) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA_Mesenchymal.pdf", useDingbats=FALSE)


#### Proneural ####
PN = read.table("Proneural - Pathways to plot.csv", sep=";", header = TRUE)
PN$N = c(1:nrow(PN))

p = ggplot(PN, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(PN$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(-max(abs(PN$Enrichment)), max(abs(PN$Enrichment))) +
  geom_vline(xintercept = 0) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA_Proneural.pdf", useDingbats=FALSE)


#### GREAT ####
setwd("D:/R/GBM-LGG/GBM Subtypes/GREAT")

#### Classical ####
CLup = read.table("GREAT Classical Up.csv", sep=";", header = TRUE, row.names=1)
CLup$N = c(1:nrow(CLup))

p = ggplot(CLup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(CLup$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(CLup$Enrichment))) +
  theme_minimal(base_size = 20)
p
ggsave(plot=p,height=NA,width=NA, filename="GREAT Classical Up.pdf", useDingbats=FALSE)

CLdown = read.table("GREAT Classical Down.csv", sep=";", header = TRUE, row.names=1)
CLdown$N = c(1:nrow(CLdown))

p = ggplot(CLdown, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(CLdown$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(CLdown$Enrichment))) +
  theme_minimal(base_size = 20)
p
ggsave(plot=p,height=NA,width=NA, filename="GREAT Classical Down.pdf", useDingbats=FALSE)

df = data.frame(Enrichment=1, ev2=c(5, 20, 40, 60), ev3=1:4)
pd = ggplot(df, aes(x=Enrichment, y=ev3, color = -5)) + 
  geom_point(size=(df$ev2*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  theme_minimal(base_size = 20)
pd
ggsave(plot=pd,height=NA,width=NA, filename="Points_for_legend.pdf", useDingbats=FALSE)



#### Mesenchymal ####
MEup = read.table("GREAT Mesenchymal Up.csv", sep=";", header = TRUE, row.names=1)
MEup$N = c(1:nrow(MEup))

p = ggplot(MEup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(MEup$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(MEup$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Mesenchymal Up.pdf", useDingbats=FALSE)


MEdown = read.table("GREAT Mesenchymal Down.csv", sep=";", header = TRUE, row.names=1)
MEdown$N = c(1:nrow(MEdown))

p = ggplot(MEdown, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(MEdown$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(MEdown$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Mesenchymal Down.pdf", useDingbats=FALSE)

#### Proneural ####
PNup = read.table("GREAT Proneural Up.csv", sep=";", header = TRUE, row.names=1)
PNup$N = c(1:nrow(PNup))

p = ggplot(PNup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(PNup$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(PNup$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Proneural Up.pdf", useDingbats=FALSE)


PNdown = read.table("GREAT Proneural Down.csv", sep=";", header = TRUE, row.names=1)
PNdown$N = c(1:nrow(PNdown))

p = ggplot(PNdown, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(PNdown$Number.of.genes*5/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(PNdown$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Proneural Down.pdf", useDingbats=FALSE)



#### GSEA SEPARATING UP/DOWN ####
#### Classical ####
CLup = CL[4:6,]
CLup$N = c(1:nrow(CLup))

p = ggplot(CLup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(CLup$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(CLup$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Classical Up.pdf", useDingbats=FALSE)


df = data.frame(Enrichment=1, ev2=c(40, 80, 120), ev3=1:3)
pd = ggplot(df, aes(x=Enrichment, y=ev3, color = -5)) + 
  geom_point(size=(df$ev2*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  theme_minimal(base_size = 20)
pd
ggsave(plot=pd,height=NA,width=NA, filename="Points_for_legend.pdf", useDingbats=FALSE)

CLdown = CL[1:3,]
CLdown$N = c(1:nrow(CLdown))

p = ggplot(CLdown, aes(x=-Enrichment, y=N, color = LogP)) + 
  geom_point(size=(CLdown$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(CLdown$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Classical Down.pdf", useDingbats=FALSE)


#### Mesenchymal ####
MEup = ME[4:6,]
MEup$N = c(1:nrow(MEup))

p = ggplot(MEup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(MEup$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(MEup$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Mesenchymal Up.pdf", useDingbats=FALSE)


MEdown = ME[1:3,]
MEdown$N = c(1:nrow(MEdown))

p = ggplot(MEdown, aes(x=-Enrichment, y=N, color = LogP)) + 
  geom_point(size=(MEdown$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(MEdown$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Mesenchymal Down.pdf", useDingbats=FALSE)

#### Proneural ####
PNup = PN[4:6,]
PNup$N = c(1:nrow(PNup))

p = ggplot(PNup, aes(x=Enrichment, y=N, color = LogP)) + 
  geom_point(size=(PNup$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(PNup$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Proneural Up.pdf", useDingbats=FALSE)


PNdown = PN[1:3,]
PNdown$N = c(1:nrow(PNdown))

p = ggplot(PNdown, aes(x=-Enrichment, y=N, color = LogP)) + 
  geom_point(size=(PNdown$GeneInGOAndHitList*3/10)) +
  labs(x = "Enrichment",
       y = "N",
       color = "LogP") +
  xlim(0, max(abs(PNdown$Enrichment))) +
  theme_minimal(base_size = 20)
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Proneural Down.pdf", useDingbats=FALSE)


#### SUPPLEMENTARY FIGURES ####
#### GSEA using 26 pathways ####

#### Classical ####
setwd("D:/R/GBM-LGG/GBM Subtypes/GSEA")
CL = read.table("Classical - 26 pathways GSEA.csv", sep=";", header = TRUE, row.names=1)
CL$N = c(1:nrow(CL))

axisLabels.x <- rownames(CL)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(CL)], y=CL$LogP, z=CL$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Classical 26 Pathways.pdf", useDingbats=FALSE)

#### Mesenchymal ####
CL = read.table("Mesenchymal - 26 pathways GSEA.csv", sep=";", header = TRUE, row.names=1)
CL$N = c(1:nrow(CL))

axisLabels.x <- rownames(CL)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(CL)], y=CL$LogP, z=CL$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Mesenchymal 26 Pathways.pdf", useDingbats=FALSE)

#### Proneural ####
CL = read.table("Proneural - 26 pathways GSEA.csv", sep=";", header = TRUE, row.names=1)
CL$N = c(1:nrow(CL))

axisLabels.x <- rownames(CL)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(CL)], y=CL$LogP, z=CL$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GSEA Proneural 26 Pathways.pdf", useDingbats=FALSE)

#### GREAT using all pathways ####
#### Classical ####
CL = read.table("Great_Classical_All_Pathways.csv", sep=";", header = TRUE, row.names=1)
CL$LogP = -log10(CL$LogP)
CL$N = c(1:nrow(CL))

axisLabels.x <- rownames(CL)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(CL)], y=CL$LogP, z=CL$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Classical All Pathways.pdf", useDingbats=FALSE)

#### Mesenchymal ####
ME = read.table("Great_Mesenchymal_All_Pathways.csv", sep=";", header = TRUE, row.names=1)
ME$LogP = -log10(ME$LogP)
ME$N = c(1:nrow(ME))

axisLabels.x <- rownames(ME)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(ME)], y=ME$LogP, z=ME$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Mesenchymal All Pathways.pdf", useDingbats=FALSE)

#### Proneural ####
PN = read.table("Great_Proneural_All_Pathways.csv", sep=";", header = TRUE, row.names=1)
PN$LogP = -log10(PN$LogP)
PN$N = c(1:nrow(PN))

axisLabels.x <- rownames(PN)
labels.wrap  <- lapply(strwrap(axisLabels.x,50,simplify=F),paste,collapse="\n") # word wrap
gg2 <- data.frame(x=LETTERS[1:nrow(PN)], y=PN$LogP, z=PN$Up.Down)


p = ggplot(gg2) +
  geom_bar(aes(x,y, fill = z), stat="identity", show.legend = TRUE)+
  guides(fill = guide_legend(title = "", title.position = "bottom")) +
  theme_minimal(base_size = 10) +
  scale_x_discrete(labels=labels.wrap)+
  scale_fill_discrete(guide="none")+
  labs(x="",y="LogP")+
  coord_flip() 
p

ggsave(plot=p,height=NA,width=NA, filename="GREAT Proneural All Pathways.pdf", useDingbats=FALSE)

