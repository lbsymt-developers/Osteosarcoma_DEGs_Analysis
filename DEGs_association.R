library(VennDiagram)
library(tidyverse)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(pathfindR)


#CARGAR LOS ARCHIVOS
cancer <- read.csv("Cancer_DEGs_Anotado.txt", fill = TRUE, header = TRUE, 
                    stringsAsFactors = F, blank.lines.skip = TRUE, row.names = 1)

saos <- read.csv("saos_ED_anotado.csv", fill = TRUE, header = TRUE, 
                stringsAsFactors = F, blank.lines.skip = TRUE, row.names = 1)

sjsa <- read.csv("sjsa_ED_anotado.csv", fill = TRUE, header = TRUE, 
                stringsAsFactors = F, blank.lines.skip = TRUE, row.names = 1)


# Venn Diagrams
cancer[cancer$refseq!="",1]
saos[saos$refseq!="",1]
sjsa[sjsa$refseq!="",1]


#Asociacion
Osteosarcoma <- list(Saos2=saos[saos$refseq!="",1], SJSA1=sjsa[sjsa$refseq!="",1], 
                     Cancer=cancer[cancer$refseq!="",1])

venn.diagram(Osteosarcoma,
             filename="Asoc_Osteosarcoma.png",
             col= "transparent", 
             fill = c("green", "blue", "red"),
             alpha=0.5, imagetype="png", 
             resolution = 300, fontfamily = "sans", 
             cat.fontfamily = "sans", height = 800, 
             width = 800, cex=1, cat.cex=1, 
             cat.col = "black", cat.pos= 1)


Asociacion1 <- calculate.overlap(list(Saos2=saos[saos$refseq!="",1], SJSA1=sjsa[sjsa$refseq!="",1], 
                                      Cancer=cancer[cancer$refseq!="",1]))

write.csv(Asociacion1$a5, file = "Asoc_all_common.csv")

#===================== IDENTIFICACION DE GENES ============================

all_cell <- as.data.frame(Asociacion1$a2)
colnames(all_cell) <- c("refseq")
all_cell2 <- merge(all_cell, sjsa, by.x = "refseq", by.y = "refseq", all.x = TRUE)

colnames(all_cell2) <- c("refseq", "logFC", "logCPM", "PValue", "FDR", "gene", "variant")
colnames(all_os) <- c("refseq", "logFC", "logCPM", "PValue", "FDR", "gene", "variant")

all_ok <- rbind(all_cell2, all_os)

all <- as.data.frame(Asociacion1$a5)
colnames(all) <- c("refseq")
all_os <- merge(all, cancer, by.x = "refseq", by.y = "refseq", all.x = TRUE)

common_ca <- as.data.frame(cbind(all_ok$gene, all_ok$logFC, all_ok$PValue))
colnames(common_ca) <- c("Gene.symbol", "logFC", "adj.P.Val")
common_ca <- common_ca[!duplicated(common_ca$Gene.symbol), ]
common_ca$logFC <- as.numeric(common_ca$logFC)
common_ca$adj.P.Val <- as.numeric(common_ca$adj.P.Val)

#KEGG
common_resultKEGG <- run_pathfindR(common_ca, gene_sets = "KEGG")
write.csv(common_resultKEGG, "KEGG_Common_OS.csv")


enrichment_chart(common_resultKEGG, top_terms = 15,
                 plot_by_cluster = FALSE,
                 num_bubbles = 3,
                 even_breaks = TRUE)

png("KEGG_common_OS.png", width=2000, height=1500, units="px", pointsize=12, res=300)



kegg_os <- read.delim("kegg_top15.txt", sep = "\t", header = TRUE)


enrichment_chart(kegg_os, top_terms = 15,
                 plot_by_cluster = FALSE,
                 num_bubbles = 3,
                 even_breaks = TRUE)

png("KEGG_common_OS_all.png", width=2000, height=1200, units="px", pointsize=12, res=300)



#GO-BP
common_resultGO-BP <- run_pathfindR(common_ca, gene_sets = "GO-BP")
write.csv(common_resultKEGG, "KEGG_Common_OS.csv")


enrichment_chart(common_resultKEGG, top_terms = 15,
                 plot_by_cluster = FALSE,
                 num_bubbles = 3,
                 even_breaks = TRUE)

png("KEGG_common_OS.png", width=2000, height=1500, units="px", pointsize=12, res=300)


#REACTOME
common_resultReact <- run_pathfindR(common_ca, gene_sets = "Reactome")
write.csv(common_resultKEGG, "KEGG_Common_OS.csv")
