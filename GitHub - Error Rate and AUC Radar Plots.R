library(fmsb)

#### AUC ####
dataAUC =as.data.frame(cbind(100, 100, 100))
dataAUC = rbind(dataAUC, c(90, 90, 90))
dataAUC = rbind(dataAUC, c(97.5, 95, 95.3))
dataAUC = rbind(dataAUC, c(90.5, 90.5, 94.8))
dataAUC = rbind(dataAUC, c(94.2, 91.4, 95.2))
colnames(dataAUC) = c("Classical", "Mesenchymal", "Proneural")
rownames(dataAUC) = c("Min", "Max", "Integrative", "Gene Expression", "DNA Methylation")

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

pdf("AUC - Radar Plot.pdf")
# plot with default options:
radarchart( dataAUC  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(90,100,2.5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(dataAUC[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()


#### ERROR RATE ####
dataER =as.data.frame(cbind(0.15, 0.15, 0.15))
dataER = rbind(dataER, c(0, 0, 0))
dataER = rbind(dataER, c((0.22-0.074), (0.22-0.106), (0.22-0.104)))
dataER = rbind(dataER, c((0.22-0.123), (0.22-0.2), (0.22-0.103)))
dataER = rbind(dataER, c((0.22-0.117), (0.22-0.163), (0.22-0.073)))
colnames(dataER) = c("Classical", "Mesenchymal", "Proneural")
rownames(dataER) = c("Min", "Max", "Integrative", "Gene Expression", "DNA Methylation")


pdf("Error Rate - All Subtypes and Methodologies in One Plot.pdf")
# plot with default options:
radarchart( dataER  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0.07,0.22,(0.15/4)), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(dataER[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
dev.off()
