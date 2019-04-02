library(tidyverse)
library(UpSetR)


currentdate <- "2019-02-15"
maxyear <- 2018


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
omim.monogenic <- subset(omim, subset=(phenomappingkey!=4&isComplex!="yes"&isComplex!="cancer"&yearDelineated<=maxyear&yearDiscovered<=maxyear))
omim.monogenic.solved <- subset(omim.monogenic, subset=(phenomappingkey==3))

omim_genes <- unique(omim.monogenic.solved$locussymbol)
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

write.csv(combined, "combined_genes.tsv", row.names=FALSE)

colSums(combined)


upset(combined, order.by="freq", nsets=8, group.by="sets")


#######################################
# trying a couple way of drawing Venn diagrams. Need to do some post-processing in Illustrator to make it pretty
#######################################


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

# qualitativecolors <- as.character(c("#332288", "#88CCEE", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"))
rainbowcolors <- as.character(c("#781C81", "#3F56A7", "#4B91C0", "#5FAA9F", "#91BD61", "#D8AF3D", "#E77C30", "#D92120"))
mycolors.red <- "#D92120"
mycolors.purple <- "#404096"
mycolors.orange <- "#E68B33"
mycolors.teal <- "#63AD99"
mycolors.blue <- "#498CC2"
mycolors.yellow <- "#BEBC48"


showSVG(outFile="Novel_known_mouse.svg", nVennObj = nVennR_mouse, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.blue), labelRegions=FALSE)
showSVG(outFile="Novel_known_human.svg", nVennObj = nVennR_human, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.orange), labelRegions=FALSE)
showSVG(outFile="Novel_known_union.svg", nVennObj = nVennR_union, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.teal), labelRegions=FALSE)
showSVG(outFile="Novel_known_intersect.svg", nVennObj = nVennR_intersect, opacity = 0.1, borderWidth = 3, setColors=c(mycolors.red,mycolors.purple,mycolors.yellow), labelRegions=FALSE)
