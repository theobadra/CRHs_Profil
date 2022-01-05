library(igraph)
library(rtracklayer)
#Exemple creation CRHs
setwd("C:\\Users\\loicm\\Documents\\DataTrek\\Mentorat")

GM12878 = read.table("GM12878.PositivePredictions.txt", header=T)

GM12878$name = paste0(GM12878$chr, ":" ,GM12878$start, ":" ,GM12878$end)

#Creation des CRHs
Enhancers.GM12878 = unique(GM12878$name)
Genes.GM12878 = unique(GM12878$TargetGene)

Enhancers.Genes.GM12878 = data.frame("Element"=rbind(matrix(Enhancers.GM12878, ncol=1), matrix(Genes.GM12878, ncol=1)))
Enhancers.Genes.GM12878$type = ifelse(Enhancers.Genes.GM12878$Element%in%Enhancers.GM12878, "Enhancer", "Gene")
Enhancers.Genes.GM12878$col = ifelse(Enhancers.Genes.GM12878$type=="Enhancer", "red", "blue")
contacts.E.G.GM12878 = GM12878[,c("name", "TargetGene")]

graph.GM12878 = graph_from_data_frame(contacts.E.G.GM12878, directed = F, vertices=Enhancers.Genes.GM12878)
components(graph.GM12878)$membership

decompose.GM12878 = decompose(graph.GM12878)
plot(decompose.GM12878[[4]], vertex.label=NA, vertex.color=V(decompose.GM12878[[4]])$col)

CTCF.narrow = import("ENCFF796WRU.bed", format = "narrowpeak")
H3K27me.narrow = import("ENCFF523KGZ.bed", format="narrowpeak")
H3K27ac.narrow = import("ENCFF411MHX.bed", format="narrowpeak")
H3K4me1.narrow = import("ENCFF507SNM.bed", format="narrowpeak")
EP300.narrow = import("ENCFF476RII.bed", format="narrowpeak")
H3K9me3.narrow = import("ENCFF713DDN.bed", format="narrowpeak")
POLR2A.narrow = import("ENCFF587YZD.bed", format="narrowpeak")
H3K36me3.narrow = import("ENCFF348OWQ.bed", format="narrowpeak")
H3K4me3.narrow = import("ENCFF656JPY.bed", format="narrowpeak")
DNAse.narrow = import("ENCSR000EMT_rep2_1_se_bwa_biorep_filtered_peaks.bed", format="narrowpeak")

GRanges.Enhancers = unique(GRanges(seqnames = GM12878$chr, ranges = IRanges(start=GM12878$start,end=GM12878$end, name=GM12878$name)))

Overlap.Enhancers.CTCF = findOverlaps(GRanges.Enhancers,CTCF.narrow)
subset.CTCF = CTCF.narrow[subjectHits(Overlap.Enhancers.CTCF)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.CTCF)), "mean_CTCF"] = aggregate(subset.CTCF$signalValue, list(queryHits(Overlap.Enhancers.CTCF)), mean, na.rm=T)$x

Overlap.Enhancers.H3K27me = findOverlaps(GRanges.Enhancers,H3K27me.narrow)
subset.H3K27me = H3K27me.narrow[subjectHits(Overlap.Enhancers.H3K27me)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K27me)), "mean_H3K27me"] = aggregate(subset.H3K27me$signalValue, list(queryHits(Overlap.Enhancers.H3K27me)), mean, na.rm=T)$x

Overlap.Enhancers.H3K27ac = findOverlaps(GRanges.Enhancers,H3K27ac.narrow)
subset.H3K27ac = H3K27ac.narrow[subjectHits(Overlap.Enhancers.H3K27ac)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K27ac)), "mean_H3K27ac"] = aggregate(subset.H3K27ac$signalValue, list(queryHits(Overlap.Enhancers.H3K27ac)), mean, na.rm=T)$x

Overlap.Enhancers.H3K4me1 = findOverlaps(GRanges.Enhancers,H3K4me1.narrow)
subset.H3K4me1 = H3K4me1.narrow[subjectHits(Overlap.Enhancers.H3K4me1)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K4me1)), "mean_H3K4me1"] = aggregate(subset.H3K4me1$signalValue, list(queryHits(Overlap.Enhancers.H3K4me1)), mean, na.rm=T)$x

Overlap.Enhancers.EP300 = findOverlaps(GRanges.Enhancers,EP300.narrow)
subset.EP300 = EP300.narrow[subjectHits(Overlap.Enhancers.EP300)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.EP300)), "mean_EP300"] = aggregate(subset.EP300$signalValue, list(queryHits(Overlap.Enhancers.EP300)), mean, na.rm=T)$x

Overlap.Enhancers.H3K9me3 = findOverlaps(GRanges.Enhancers,H3K9me3.narrow)
subset.H3K9me3 = H3K9me3.narrow[subjectHits(Overlap.Enhancers.H3K9me3)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K9me3)), "mean_H3K9me3"] = aggregate(subset.H3K9me3$signalValue, list(queryHits(Overlap.Enhancers.H3K9me3)), mean, na.rm=T)$x

Overlap.Enhancers.POLR2A = findOverlaps(GRanges.Enhancers,POLR2A.narrow)
subset.POLR2A = POLR2A.narrow[subjectHits(Overlap.Enhancers.POLR2A)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.POLR2A)), "mean_POLR2A"] = aggregate(subset.POLR2A$signalValue, list(queryHits(Overlap.Enhancers.POLR2A)), mean, na.rm=T)$x

Overlap.Enhancers.H3K36me3 = findOverlaps(GRanges.Enhancers,H3K36me3.narrow)
subset.H3K36me3 = H3K36me3.narrow[subjectHits(Overlap.Enhancers.H3K36me3)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K36me3)), "mean_H3K36me3"] = aggregate(subset.H3K36me3$signalValue, list(queryHits(Overlap.Enhancers.H3K36me3)), mean, na.rm=T)$x

Overlap.Enhancers.H3K4me3 = findOverlaps(GRanges.Enhancers,H3K4me3.narrow)
subset.H3K4me3 = H3K4me3.narrow[subjectHits(Overlap.Enhancers.H3K4me3)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.H3K4me3)), "mean_H3K4me3"] = aggregate(subset.H3K4me3$signalValue, list(queryHits(Overlap.Enhancers.H3K4me3)), mean, na.rm=T)$x

Overlap.Enhancers.DNAse = findOverlaps(GRanges.Enhancers,DNAse.narrow)
subset.DNAse = DNAse.narrow[subjectHits(Overlap.Enhancers.DNAse)]
mcols(GRanges.Enhancers)[unique(queryHits(Overlap.Enhancers.DNAse)), "mean_DNAse"] = aggregate(subset.DNAse$signalValue, list(queryHits(Overlap.Enhancers.DNAse)), mean, na.rm=T)$x

library(reshape2)
library(FactoMineR)
library(dplyr)

epi.signal = data.frame(mcols(GRanges.Enhancers))
epi.signal[is.na(epi.signal)] = 0
epi.signal$name = rownames(epi.signal)

df.membership = melt(components(graph.GM12878)$membership)
df.membership$name = rownames(df.membership)

epi.signal.membership = merge(epi.signal, df.membership, by="name")
epi.signal.membership$membership = as.factor(epi.signal.membership$value)

epi.signal.agg.membership = data.frame(epi.signal.membership %>% group_by(membership)%>%summarise_at(vars("mean_H3K4me3",
                                                                "mean_DNAse", "mean_CTCF", "mean_H3K27me",
                                                                "mean_H3K27ac", "mean_H3K4me1", "mean_EP300", "mean_H3K9me3",
                                                                "mean_POLR2A", "mean_H3K36me3"), mean))

res.pca = PCA(epi.signal.membership[,2:10], graph=T, scale.unit = T)
plot.PCA(res.pca, axes=c(1, 2), choix="ind", habillage="none", label = "none")

