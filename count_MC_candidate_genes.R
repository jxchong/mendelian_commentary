library(tidyverse)

currentdate <- "2019-02-15"


# read in data
raw_hgncgenes <- read.table("HGNC_genes.tsv", head=FALSE, stringsAsFactors=FALSE)
allgenes <- raw_hgncgenes[,1]

raw_NMDesc <- read.table("NMDescape.txt", head=TRUE, stringsAsFactors=FALSE)
raw_ccrs_autosomes <- read.table("ccrs.autosomes.90orhigher.v2.20180420.bed.gz", head=TRUE, comment.char="", stringsAsFactors=FALSE)
raw_ccrs_xchr <- read.table("ccrs.xchrom.90orhigher.v2.20180420.bed.gz", head=TRUE, comment.char="", stringsAsFactors=FALSE)

raw_gnomad_constraint <- read.table("gnomad_constraint.txt.bgz", head=TRUE, stringsAsFactors=FALSE)
raw_gnomad_canonical <- raw_gnomad_constraint[raw_gnomad_constraint$canonical=="true", c("gene", "oe_mis_upper", "oe_lof_upper", "mis_z", "pLI")]
raw_gnomad_canonical$mis_rank <- NA
raw_gnomad_canonical$lof_rank <- NA
raw_gnomad_canonical$mis_rank[order(raw_gnomad_canonical$oe_mis_upper)] <- 1:nrow(raw_gnomad_canonical)
raw_gnomad_canonical$lof_rank[order(raw_gnomad_canonical$oe_lof_upper)] <- 1:nrow(raw_gnomad_canonical)
raw_gnomad_canonical_mis <- raw_gnomad_canonical %>% mutate(mis_ile = ntile(mis_rank, 10))
raw_gnomad_canonical_mis_lof <- raw_gnomad_canonical_mis %>% mutate(lof_ile = ntile(lof_rank, 10))

omim <- read.table(paste0(currentdate,".combinedOMIM.mentionsNGS.year.inheritance.txt"), head=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", strip.white=TRUE)
omim.monogenic <- subset(omim, subset=(phenomappingkey!=4&isComplex!="yes"&isComplex!="cancer"))
omim.monogenic.solved <- subset(omim.monogenic, subset=(phenomappingkey==3))
omim_genes <- unique(omim.monogenic$locussymbol)
potentialMC_genes <- unique(allgenes[!(allgenes %in% omim_genes)])

ccr_genes <- unique(c(raw_ccrs_autosomes$gene, raw_ccrs_xchr$gene))
gnomad_mis_genes <- unique(raw_gnomad_canonical_mis_lof[raw_gnomad_canonical_mis_lof$mis_ile<=2,"gene"])
gnomad_lof_genes <- unique(raw_gnomad_canonical_mis_lof[raw_gnomad_canonical_mis_lof$lof_ile<=4,"gene"])
gnomad_pli_genes <- unique(raw_gnomad_canonical_mis_lof[raw_gnomad_canonical_mis_lof$pLI>=0.9,"gene"])


# read in mouse data
raw_mouse <- read.table("HMD_HumanPhenotype.rpt.mortality_HPOcount.txt", head=TRUE, stringsAsFactors=FALSE, sep="\t")
mouse_lethal_genes <- unique(raw_mouse[raw_mouse$MortalityBool==TRUE, "Human.Marker.Symbol"])
mouse_nonlethal_genes <- unique(raw_mouse[raw_mouse$HPOcount>(raw_mouse$NormalBool+raw_mouse$NoAnalysisBool+raw_mouse$MortalityBool), "Human.Marker.Symbol"])
mouse_pheno_genes <- unique(c(mouse_nonlethal_genes, mouse_lethal_genes))
mouse_notstudied_genes <- unique(raw_mouse[raw_mouse$HPOcount==0|raw_mouse$NoAnalysisBool==TRUE, "Human.Marker.Symbol"])
mouse_normal_genes <- unique(raw_mouse[raw_mouse$HPOcount==1&raw_mouse$NormalBool==TRUE, "Human.Marker.Symbol"])
mouse_human_genes <- unique(raw_mouse$Human.Marker.Symbol)

combined <- data.frame(allgenes,row.names=1,stringsAsFactors=FALSE)
combined$ccr <- 0
combined$NMDesc <- 0
combined$gnomad_mis <- 0
combined$gnomad_lof <- 0
combined$any_predictor <- 0
combined$OMIMknown <- 0
combined$potentialMC <- 0
combined$mouse_lethal <- 0
combined$mouse_nonlethal <- 0
combined$mouse_pheno <- 0
combined$humanmouse_union <- 0
combined$humanmouse_intersect <- 0
combined[row.names(combined) %in% raw_NMDesc$gene, "NMDesc"] <- 1
combined[row.names(combined) %in% gnomad_mis_genes, "gnomad_mis"] <- 1
combined[row.names(combined) %in% gnomad_lof_genes, "gnomad_lof"] <- 1
combined$any_predictor <- as.integer((combined$ccr + combined$NMDesc + combined$gnomad_mis + combined$gnomad_lof) >= 1)
combined[row.names(combined) %in% omim_genes, "OMIMknown"] <- 1
combined[row.names(combined) %in% potentialMC_genes, "potentialMC"] <- 1
combined[row.names(combined) %in% mouse_lethal_genes, "mouse_lethal"] <- 1
combined[row.names(combined) %in% mouse_nonlethal_genes, "mouse_nonlethal"] <- 1
combined[row.names(combined) %in% mouse_pheno_genes, "mouse_pheno"] <- 1
combined$humanmouse_union <- as.integer((combined$any_predictor + combined$mouse_pheno) >= 1)
combined$humanmouse_intersect <- as.integer((combined$any_predictor+combined$mouse_pheno) == 2)

write.csv(combined, "combined_genes.tsv", row.names=TRUE)

colSums(combined)
# > colSums(combined)
#                  ccr               NMDesc           gnomad_mis           gnomad_lof        any_predictor            OMIMknown          potentialMC
#                    0                 1890                 3743                 7540                 9596                 3519                15675
#         mouse_lethal      mouse_nonlethal          mouse_pheno     humanmouse_union humanmouse_intersect
#                 5497                10199                10487                13737                 6346


# some data exploration
# what are the characteristics of the genes with no human or mouse support? genes with no human support?
# answer: not surprisingly, it's enriched for recessives

solvednohumanmouse_union <- subset(omim.monogenic, locussymbol%in%allgenes[combined$humanmouse_union==FALSE&combined$OMIMknown==TRUE])
solvednohumanmouse_union.annot <- merge(solvednohumanmouse_union, raw_gnomad_canonical_mis_lof, by.x="locussymbol", by.y="gene")

solvednohumanmouse_union.inheritDN <- subset(solvednohumanmouse_union.annot, subset=(official_inheritDN==1|inheritDN==1)&(official_inheritAD==1|inheritAD==1|official_inheritX==1|inheritX==1)&mentionsLinkage==FALSE)
solvednohumanmouse_union.inheritAD <- subset(solvednohumanmouse_union.annot, subset=(official_inheritAD==1|inheritAD==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE)&(official_inheritAR==0|inheritAR==0))
solvednohumanmouse_union.inheritAR <- subset(solvednohumanmouse_union.annot, subset=(official_inheritAR==1|inheritAR==1)&(official_inheritAD==0|inheritAD==0))
solvednohumanmouse_union.inheritX <- subset(solvednohumanmouse_union.annot, subset=(official_inheritX==1|inheritX==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE))

layout(matrix(1:8, ncol=2, byrow=TRUE))
hist(solvednohumanmouse_union.inheritDN$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritDN$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritAD$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritAD$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritAR$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritAR$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritX$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohumanmouse_union.inheritX$lof_rank, xlim=c(1,20000), ylim=c(0,200))


# store omim info for known genes with no human constraint support
solvednohuman <- subset(omim.monogenic, locussymbol%in%allgenes[combined$any_predictor==FALSE&combined$OMIMknown==TRUE])
solvednohuman.annot <- merge(solvednohuman, raw_gnomad_canonical_mis_lof, by.x="locussymbol", by.y="gene")
write.table(solvednohuman.annot, file="knowngenes_nohumansupport.tsv", quote=FALSE, sep="\t", row.names=FALSE)

solvednohuman.inheritDN <- subset(solvednohuman.annot, subset=(official_inheritDN==1|inheritDN==1)&(official_inheritAD==1|inheritAD==1|official_inheritX==1|inheritX==1)&mentionsLinkage==FALSE)
solvednohuman.inheritAD <- subset(solvednohuman.annot, subset=(official_inheritAD==1|inheritAD==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE)&(official_inheritAR==0|inheritAR==0))
solvednohuman.inheritAR <- subset(solvednohuman.annot, subset=(official_inheritAR==1|inheritAR==1)&(official_inheritAD==0|inheritAD==0))
solvednohuman.inheritX <- subset(solvednohuman.annot, subset=(official_inheritX==1|inheritX==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE))

layout(matrix(1:8, ncol=2, byrow=TRUE))
hist(solvednohuman.inheritDN$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritDN$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritAD$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritAD$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritAR$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritAR$lof_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritX$mis_rank, xlim=c(1,20000), ylim=c(0,200))
hist(solvednohuman.inheritX$lof_rank, xlim=c(1,20000), ylim=c(0,200))





#######################################
#  various plotting methods
#######################################

# library(UpSetR) (for interactive exploratory visualization)
upset(combined, order.by="freq", nsets=8, group.by="sets")


# qualitativecolors <- as.character(c("#332288", "#88CCEE", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"))
rainbowcolors <- as.character(c("#781C81", "#3F56A7", "#4B91C0", "#5FAA9F", "#91BD61", "#D8AF3D", "#E77C30", "#D92120"))
mycolors.red <- "#D92120"
mycolors.purple <- "#404096"
mycolors.orange <- "#E68B33"
mycolors.teal <- "#63AD99"
mycolors.blue <- "#498CC2"
mycolors.yellow <- "#BEBC48"




library(eulerr)

# error_plot(v)

pdf("Potential_novel_Mendelian_genes_and_known.pdf", width=6, height=6)
plot(euler(
    c("MouseGenesPheno"=10487, "HumanConstrainedGenes"=13190, "PotentialMendelianGenes"=15675, "KnownMendelianGenes"=3519,
    "MouseGenesPheno&HumanConstrainedGenes"=8302,
    "MouseGenesPheno&PotentialMendelianGenes"=7501,
    "MouseGenesPheno&KnownMendelianGenes"=2986,
    "HumanConstrainedGenes&PotentialMendelianGenes"=10276,
    "HumanConstrainedGenes&KnownMendelianGenes"=2914,
    "PotentialMendelianGenes&KnownMendelianGenes"=0,
    "MouseGenesPheno&HumanConstrainedGenes&PotentialMendelianGenes"=5791,
    "MouseGenesPheno&HumanConstrainedGenes&KnownMendelianGenes"=2511,
    "HumanConstrainedGenes&PotentialMendelianGenes&KnownMendelianGenes"=0,
    "MouseGenesPheno&HumanConstrainedGenes&PotentialMendelianGenes&KnownMendelianGenes"=0
    ),
    input = "union",
    shape = "ellipse"
    ),
    quantities = TRUE,
    legend=TRUE
)
dev.off()


pdf("Potential_novel_Mendelian_genes_and_known.circle.pdf", width=6, height=6)
plot(euler(
    c("MouseGenesPheno"=10487, "HumanConstrainedGenes"=13190, "PotentialMendelianGenes"=15675, "KnownMendelianGenes"=3519,
    "MouseGenesPheno&HumanConstrainedGenes"=8302,
    "MouseGenesPheno&PotentialMendelianGenes"=7501,
    "MouseGenesPheno&KnownMendelianGenes"=2986,
    "HumanConstrainedGenes&PotentialMendelianGenes"=10276,
    "HumanConstrainedGenes&KnownMendelianGenes"=2914,
    "PotentialMendelianGenes&KnownMendelianGenes"=0,
    "MouseGenesPheno&HumanConstrainedGenes&PotentialMendelianGenes"=5791,
    "MouseGenesPheno&HumanConstrainedGenes&KnownMendelianGenes"=2511,
    "HumanConstrainedGenes&PotentialMendelianGenes&KnownMendelianGenes"=0,
    "MouseGenesPheno&HumanConstrainedGenes&PotentialMendelianGenes&KnownMendelianGenes"=0
    ),
    input = "union"
    ),
    quantities = TRUE,
    lty=1:4,
    legend=TRUE
)
dev.off()

# > colSums(combined)
#                  ccr               NMDesc           gnomad_mis           gnomad_lof        any_predictor            OMIMknown          potentialMC
#                    0                 1890                 3743                 7540                 9596                 3519                15675
#         mouse_lethal      mouse_nonlethal          mouse_pheno     humanmouse_union humanmouse_intersect
#                 5497                10199                10487                13737                 6346


euler.union <- euler(
    c("Candidates"=0, "KnownMendelianGenes"=249, "PotentialMendelianGenes"=5208,
    "Candidates&KnownMendelianGenes"=3270,
    "Candidates&PotentialMendelianGenes"=10467,
    "Candidates&KnownMendelianGenes&PotentialMendelianGenes"=0
    ),
    input = "disjoint",
    shape = "circle"
)
euler.intersect <- euler(
    c("CandidatesIntersect"=0, "KnownMendelianGenes"=1623, "PotentialMendelianGenes"=11225,
    "CandidatesIntersect&KnownMendelianGenes"=1896,
    "CandidatesIntersect&PotentialMendelianGenes"=4450,
    "CandidatesIntersect&KnownMendelianGenes&PotentialMendelianGenes"=0
    ),
    input = "disjoint",
    shape = "circle"
)

pdf("Potential_novel_Mendelian_genes_and_known_vs_unioncandidates.circle.pdf", width=8.5, height=11)
plot(euler.union,
    quantities = TRUE
)
plot(euler.intersect,
    quantities = TRUE
)
dev.off()



test<- venn(c("CandidatesIntersect"=0, "KnownMendelianGenes"=1623, "PotentialMendelianGenes"=11225,
"CandidatesIntersect&KnownMendelianGenes"=1896,
"CandidatesIntersect&PotentialMendelianGenes"=4450,
"CandidatesIntersect&KnownMendelianGenes&PotentialMendelianGenes"=0))

#########################################

library(nVennR)

nVennR_union <- plotVenn(list(
    "Known gene"=row.names(subset(combined, OMIMknown==1)),
    "Potential novel gene"=row.names(subset(combined, potentialMC==1)),
    "Mouse or human support"=row.names(subset(combined, humanmouse_union==1))
))
nVennR_intersect <- plotVenn(list(
    "Known gene"=row.names(subset(combined, OMIMknown==1)),
    "Potential novel gene"=row.names(subset(combined, potentialMC==1)),
    "Mouse and human support"=row.names(subset(combined, humanmouse_intersect==1))
))
nVennR_mouse <- plotVenn(list(
    "Known gene"=row.names(subset(combined, OMIMknown==1)),
    "Potential novel gene"=row.names(subset(combined, potentialMC==1)),
    "Mouse support"=row.names(subset(combined, mouse_pheno==1))
))
nVennR_human <- plotVenn(list(
    "Known gene"=row.names(subset(combined, OMIMknown==1)),
    "Potential novel gene"=row.names(subset(combined, potentialMC==1)),
    "Human support"=row.names(subset(combined, any_predictor==1))
))
# listVennRegions(nVennR_union)


showSVG(outFile="Novel_known_mouse.svg", nVennObj = nVennR_mouse, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.blue), labelRegions=FALSE)
showSVG(outFile="Novel_known_human.svg", nVennObj = nVennR_human, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.orange), labelRegions=FALSE)
showSVG(outFile="Novel_known_union.svg", nVennObj = nVennR_union, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.teal), labelRegions=FALSE)
showSVG(outFile="Novel_known_intersect.svg", nVennObj = nVennR_intersect, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.yellow), labelRegions=FALSE)



### try stacked bar graph plotting


humanonly <- c("Known genes for MCs supported by human data", "Known genes for MCs not supported by human data", "Novel genes for MCs not supported by human data", "Novel genes for MCs supported by human data (priority candidate genes)")
mouseonly <- c("Known genes for MCs supported by mouse data", "Known genes for MCs not supported by mouse data", "Novel genes for MCs not supported by mouse data",  "Novel genes for MCs supported by mouse data (priority candidate genes)")
humanandmouse <- c("Known genes for MCs supported by both human and mouse data", "Known genes for MCs not supported by both", "Novel genes for MCs not supported by both", "Novel genes for MCs supported by both (priority candidate genes)")
humanormouse <- c("Known genes for MCs supported by either human or mouse data", "Known genes for MCs not supported by either", "Novel genes for MCs not supported by either", "Novel genes for MCs supported by either (priority candidate genes)")

humanormousedata <- c(249, 3270, 5208, 10467)
humanandmousedata <- c(1623, 1896, 11225, 4450)
mouseonlydata <- c(533, 2986, 8174, 7501)
humanonlydata <- c(1393, 2180, 8259, 7416)


genedata <- t(matrix(
    c(humanormousedata, humanandmousedata, mouseonlydata, humanonlydata),
    nrow=4, ncol=4, byrow=TRUE
))



# Find the top y position of each block
barmids <- apply(genedata, 2, cumsum)
# Move it downwards half the size of each block
barmids <- barmids - genedata/2
barmids <- t(barmids)



pdf("nMCs_barplot_original.pdf", width=6, height=3.3, pointsize=12)
par(mar=c(4,1,3,1.5))
bar.y <- barplot(genedata, horiz=TRUE,
    xlim=c(0,20000), xlab="Number of genes",
    col=c("#4277BD", "#529DB7", "#E39C37", "#E76D2E"), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, cex.lab=0.8
)

labelcolors <- c("white", "white", "black", "black")
text(barmids[4,], bar.y[4], lab=humanonlydata, col=labelcolors, cex=0.9)
text(barmids[3,], bar.y[3], lab=mouseonlydata, col=labelcolors, cex=0.9)
text(barmids[2,], bar.y[2], lab=humanandmousedata, col=labelcolors, cex=0.9)
text(barmids[1,], bar.y[1], lab=humanormousedata, col=labelcolors, cex=0.9)
dev.off()




text(x=bars.x, y=-20, labels=years.important, cex=0.8, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.8, line=2)
