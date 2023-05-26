library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyranges)

options(ucscChromosomeNames=FALSE)
##chr7:6,402,183-6,404,081
CHR="chr7"
START=6402183
END=6404081

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

grtrack <-  GeneRegionTrack(txdb, chromosome=CHR, start=START, end=END, symbol="Rac1",
                            showId=TRUE, name="Rac1", geneSymbol=TRUE)
tr_genome <- GenomeAxisTrack(v=10)

gr_pery <- read_bigwig("../hcl/data/bigwig/celltypes/12-erythroid-progenitor-cell-rp-high.negative.bw") %>%
    filter(seqnames == "7") %>%
    `seqlevelsStyle<-`("UCSC")

gr_endo <- read_bigwig("../hcl/data/bigwig/celltypes/29-endothelial-cell.negative.bw") %>%
    filter(seqnames == "7") %>%
    `seqlevelsStyle<-`("UCSC")

## HCL Data
tr_hcl_pery <- DataTrack(gr_pery, name="Erythroid\nProgenitor", genome="hg38",
                         type="polygon", ylim=c(0,175),col='black', area="green")
tr_hcl_endo <- DataTrack(gr_endo,
                         name="Endothelial\nCell", genome="hg38",
                         type="polygon",
                         ylim=c(0,175))

pdf("img/sc-studies/fuckoff.pdf", width=5, height=5)
plotTracks(list(tr_hcl_pery, tr_hcl_endo, grtrack, tr_genome),
           sizes=c(4,4,4,2), from=START, to=END, fontcolor='black')
dev.off()


