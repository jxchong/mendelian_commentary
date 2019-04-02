reviewred <- "#E72F1C"
reviewblue <- "#05ADEE"

phenoorange <- "#EE7733"
phenoteal <- "#009988"

cmgcolor <- "#E39C37"

# colors for venn diagrams
morecolors <- as.character(c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677", "#882255", "#AA4499"))
mycolors <- as.character(c("#D92120", "#E39C37", "#7DB874", "#529DB7"))




args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	stop("Need datestamp argument\n", call.=FALSE)
} else {
	currentdate <- paste0("", args[1])
}


# currentdate <- "2019-02-15"
maxyear <- 2017



reviewred <- "#E72F1C"
reviewblue <- "#05ADEE"
phenoorange <- "#EE7733"
phenoteal <- "#009988"



omim <- read.table(paste0(currentdate,".combinedOMIM.mentionsNGS.year.inheritance.txt"), head=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char="", quote="", strip.white=TRUE)


omim.monogenic <- subset(omim, subset=(phenomappingkey!=4&isComplex!="yes"&isComplex!="cancer"&yearDelineated<=maxyear&yearDiscovered<=maxyear))
omim.monogenic.solved <- subset(omim.monogenic, subset=(phenomappingkey==3))


geneknown <- subset(omim.monogenic, subset=yearDiscovered<=yearDelineated&phenomappingkey==3)
delineationonly <- subset(omim.monogenic, subset=((yearDiscovered>yearDelineated)|phenomappingkey<=2))

allyears <- seq(1900,maxyear,by=1)
years <- as.character(1986:maxyear)
years.important <- rep("", (maxyear-1986+1))
years.important[which(years==1986)] <- 1986
years.important[which(years==1990)] <- 1990
years.important[which(years==1995)] <- 1995
years.important[which(years==2000)] <- 2000
years.important[which(years==2005)] <- 2005
years.important[which(years==2010)] <- 2010
years.important[which(years==2015)] <- 2015


allyears <- seq(1900,maxyear,1)
allyears.important <- allyears
allyears.important[!(allyears %in% seq(1900,maxyear,10))] <- ""



# compare delineation
geneknownvsdelineation <- matrix(rbind(table(subset(delineationonly, select=yearDelineated))[years], table(subset(geneknown, select=yearDelineated))[years]), nrow=2)
colnames(geneknownvsdelineation) <- years
geneknownvsdelineation[is.na(geneknownvsdelineation)] <- 0

geneknownvsdelineation.allyears <- matrix(rbind(table(subset(delineationonly, select=yearDelineated))[as.character(allyears)], table(subset(geneknown, select=yearDelineated))[as.character(allyears)]), nrow=2)
colnames(geneknownvsdelineation.allyears) <- allyears
geneknownvsdelineation.allyears[is.na(geneknownvsdelineation.allyears)] <- 0


# 1986-1997 inclusive; 1998-2009; 2010-maxyear*
summary(lm(colSums(geneknownvsdelineation[,1:12]) ~ as.numeric(colnames(geneknownvsdelineation)[1:12]) ))
# 1998-2009
summary(lm(colSums(geneknownvsdelineation[,13:24]) ~ as.numeric(colnames(geneknownvsdelineation)[13:24]) ))
# 2010-15
summary(lm(colSums(geneknownvsdelineation[,25:30]) ~ as.numeric(colnames(geneknownvsdelineation)[25:30]) ))


geneknownvsdelineation.cumsum <- cumsum(colSums(geneknownvsdelineation))
geneknownvsdelineation.cumsum.allyears <- cumsum(colSums(geneknownvsdelineation.allyears))



current.NGS <- subset(omim.monogenic, select=NGSyear, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&mentionsNGS=="yes"&grepl("'1/", NGSparagraph))
current.noNGS <- subset(omim.monogenic, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))&yearDiscovered>=2010)
current.noNGS.old <- table(subset(omim.monogenic, select=yearDiscovered, subset=phenomappingkey==3&isComplex!="yes"&isComplex!="cancer"&!(mentionsNGS=="yes"&grepl("'1/", NGSparagraph))))

# show all years to compare NGS vs no NGS
current.NGSvsnoNGS.allyears <- as.matrix(rbind(current.noNGS.old, table(current.NGS)[names(current.noNGS.old)]))
current.NGSvsnoNGS.allyears[is.na(current.NGSvsnoNGS.allyears)] <- 0
current.NGSvsnoNGS.allyears <- current.NGSvsnoNGS.allyears[,colnames(current.NGSvsnoNGS.allyears)>=1956]
emptyyears <- matrix(0, nrow=2, ncol=length(1900:1955))
colnames(emptyyears) <- seq(1900,1955,by=1)
current.NGSvsnoNGS.allyears <- cbind(emptyyears, current.NGSvsnoNGS.allyears)

current.NGSvsnoNGS.since1986 <- current.NGSvsnoNGS.allyears[,colnames(current.NGSvsnoNGS.allyears)>=1986&colnames(current.NGSvsnoNGS.allyears)<=maxyear]
current.NGSvsnoNGS.NGSyears <- current.NGSvsnoNGS.since1986[,current.NGSvsnoNGS.since1986[2,]!=0]

# 1986-1997 inclusive; 1998-2009; 2010-2015/8
rbind(1:ncol(current.NGSvsnoNGS.since1986), current.NGSvsnoNGS.since1986)
summary(lm(colSums(current.NGSvsnoNGS.since1986[,1:12]) ~ as.numeric(colnames(current.NGSvsnoNGS.since1986)[1:12]) ))
summary(lm(colSums(current.NGSvsnoNGS.since1986[,13:24]) ~ as.numeric(colnames(current.NGSvsnoNGS.since1986)[13:24]) ))
summary(lm(colSums(current.NGSvsnoNGS.since1986[,25:30]) ~ as.numeric(colnames(current.NGSvsnoNGS.since1986)[25:30]) ))

current.NGSvsnoNGS.since1986.cumsum <- cumsum(colSums(current.NGSvsnoNGS.since1986))

# total discoveries by NGS vs not
rowSums(current.NGSvsnoNGS.since1986)






modelnames <- c("de novo", "autosomal dominant", "autosomal recessive", "X-linked")
officialmodelnames <- c("de novo", "autosomal dominant", "autosomal recessive", "X-linked")

solved.inheritDN <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&(official_inheritDN==1|inheritDN==1)&(official_inheritAD==1|inheritAD==1|official_inheritX==1|inheritX==1)&mentionsLinkage==FALSE, select=yearDiscovered))
solved.inheritAD <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&(official_inheritAD==1|inheritAD==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE)&(official_inheritAR==0|inheritAR==0), select=yearDiscovered))
solved.inheritAR <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&(official_inheritAR==1|inheritAR==1)&(official_inheritAD==0|inheritAD==0), select=yearDiscovered))
solved.inheritX <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&(official_inheritX==1|inheritX==1)&((official_inheritDN==0&inheritDN==0)|mentionsLinkage==TRUE), select=yearDiscovered))

solved.allmodels <- matrix(0, nrow=4, ncol=length(years))
colnames(solved.allmodels) <- years
rownames(solved.allmodels) <- modelnames
solved.allmodels[1,names(solved.inheritDN)] <- solved.inheritDN
solved.allmodels[2,names(solved.inheritAD)] <- solved.inheritAD
solved.allmodels[3,names(solved.inheritAR)] <- solved.inheritAR
solved.allmodels[4,names(solved.inheritX)] <- solved.inheritX


solved.official_inheritDN <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&official_inheritDN==1&(official_inheritAD==1|official_inheritX==1), select=yearDiscovered))
solved.official_inheritAD <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&official_inheritAD==1&official_inheritDN==0&official_inheritAR==0, select=yearDiscovered))
solved.official_inheritAR <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&official_inheritAR==1&official_inheritAD==0, select=yearDiscovered))
solved.official_inheritX <- table(subset(omim.monogenic.solved, subset=yearDiscovered<=maxyear&yearDiscovered>=1986&official_inheritX==1&official_inheritDN==0, select=yearDiscovered))

solved.officialallmodels <- matrix(0, nrow=4, ncol=length(years))
colnames(solved.officialallmodels) <- years
rownames(solved.officialallmodels) <- officialmodelnames
solved.officialallmodels[1,names(solved.official_inheritDN)] <- solved.official_inheritDN
solved.officialallmodels[2,names(solved.official_inheritAD)] <- solved.official_inheritAD
solved.officialallmodels[3,names(solved.official_inheritAR)] <- solved.official_inheritAR
solved.officialallmodels[4,names(solved.official_inheritX)] <- solved.official_inheritX


#  calculate de novo proportion each year
DNproportion <- mapply(`/`, solved.allmodels["de novo",],  colSums(solved.allmodels))



pdf(paste0(currentdate,".delineatedsyndromes_per_year.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(2.4,4,1,0))
bars.x <- barplot(geneknownvsdelineation, beside=FALSE, ylab="# new phenotypes reported", ylim=c(0,300), col=c(phenoorange, phenoteal), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=years.important, cex=0.7, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.7, line=1.4)
text(x=bars.x, y=colSums(geneknownvsdelineation), labels=colSums(geneknownvsdelineation), pos=3, xpd=TRUE, cex=0.5)
text(x=bars.x[geneknownvsdelineation[2,]>=10], y=geneknownvsdelineation[1,geneknownvsdelineation[2,]>=10]+0.5*geneknownvsdelineation[2,geneknownvsdelineation[2,]>=10], labels=geneknownvsdelineation[2,geneknownvsdelineation[2,]>=10], xpd=TRUE, cex=0.4, col="white")
legend("topleft", fill=c(phenoorange, phenoteal), legend=c("delineation only", "delineation and gene discovery"), border=NA, bty="n", cex=0.8)
dev.off()

# delineations older (all years)
pdf(paste0(currentdate,".delineatedsyndromes_per_year-allyears.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(3,4,1,0))
bars.x <- barplot(geneknownvsdelineation.allyears, beside=FALSE, ylab="# new phenotypes reported", ylim=c(0,300), col=c(phenoorange, phenoteal), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=allyears.important, cex=0.7, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.8, line=2)
# text(x=bars.x, y=colSums(geneknownvsdelineation.allyears), labels=colSums(geneknownvsdelineation.allyears), pos=3, xpd=TRUE, cex=0.3)
# text(x=bars.x[geneknownvsdelineation.allyears[2,]>=10], y=geneknownvsdelineation.allyears[1,geneknownvsdelineation.allyears[2,]>=10]+0.5*geneknownvsdelineation.allyears[2,geneknownvsdelineation.allyears[2,]>=10], labels=geneknownvsdelineation.allyears[2,geneknownvsdelineation.allyears[2,]>=10], xpd=TRUE, cex=0.5, col="black")
legend("topleft", fill=c(phenoorange, phenoteal), legend=c("delineation only", "delineation and gene discovery"), border=NA, bty="n", cex=0.8)
dev.off()

pdf(paste0(currentdate,".delineatedsyndromes_per_year-cumulative-1986on.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(2.8,4,1,0))
plot(1:length(geneknownvsdelineation.cumsum), geneknownvsdelineation.cumsum, ylim=c(0,5000), type="o", pch=19, cex=0.65, cex.axis=0.7, cex.lab=0.8, xaxt="n", xlab="", ylab="Cumulative new syndromes reported", bty="n")
mtext("Year", side=1, cex=0.7, line=1.9)
text(x=1:length(geneknownvsdelineation.cumsum), y=-700, labels=years.important, cex=0.7, xpd=TRUE, srt=90)
# axis(side=1, tick=TRUE, labels=FALSE, line=0.5)
points(1:length(geneknownvsdelineation.cumsum), cumsum(geneknownvsdelineation[1,]), type="p", col=phenoorange, pch=19, cex=0.65)
lines(1:length(geneknownvsdelineation.cumsum), cumsum(geneknownvsdelineation[1,]), lty="solid", col=phenoorange, pch=19)
points(1:length(geneknownvsdelineation.cumsum), cumsum(geneknownvsdelineation[2,]), type="p", col=phenoteal, pch=19, cex=0.65)
lines(1:length(geneknownvsdelineation.cumsum), cumsum(geneknownvsdelineation[2,]), lty="solid", col=phenoteal, pch=19)
legend("topleft", fill=c("black", phenoorange, phenoteal), legend=c("total new syndromes", "delineation only", "delineation and gene discovery"), border=NA, bty="n", cex=0.8)
dev.off()

pdf(paste0(currentdate,".delineatedsyndromes_per_year-cumulative-1900on.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(3,4,1,0))
plot(seq(1900,maxyear,by=1), geneknownvsdelineation.cumsum.allyears, type="o", pch=19, cex.axis=0.7, cex=0.5, cex.lab=0.8, xaxt="n", xlab="", ylab="Cumulative new syndromes reported", bty="n")
mtext("Year", side=1, cex=1, line=3)
text(seq(1900,maxyear,by=5), y=-900, labels=seq(1900,maxyear,by=5), cex=0.8, srt=90, xpd=TRUE)
points(seq(1900,maxyear,by=1), cumsum(geneknownvsdelineation.allyears[1,]), type="o", col=phenoorange, pch=19, cex=0.65)
points(seq(1900,maxyear,by=1), cumsum(geneknownvsdelineation.allyears[2,]), type="o", col=phenoteal, pch=19, cex=0.65)
legend("topleft", fill=c("black", phenoorange, phenoteal), legend=c("total new syndromes", "delineation only", "delineation and gene discovery"), border=NA, bty="n", cex=0.8)
dev.off()





pdf(paste0(currentdate,".NGSvsNoNGS_discoveries_per_year.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(2.4,4,1,0))
bars.x <- barplot(current.NGSvsnoNGS.since1986, beside=FALSE, ylab="# of gene discoveries by method", ylim=c(0,300), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=years.important, cex=0.7, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.7, line=1.4)
text(x=bars.x, y=colSums(current.NGSvsnoNGS.since1986), labels=colSums(current.NGSvsnoNGS.since1986), pos=3, xpd=TRUE, cex=0.5)
text(x=bars.x[current.NGSvsnoNGS.since1986[2,]!=0], y=current.NGSvsnoNGS.NGSyears[1,]+0.5*current.NGSvsnoNGS.NGSyears[2,], labels=current.NGSvsnoNGS.since1986[2,current.NGSvsnoNGS.since1986[2,]!=0], xpd=TRUE, cex=0.4, col="white")
legend("topleft", fill=c(reviewred, reviewblue), legend=c("conventional", "NGS"), border=NA, bty="n", cex=0.8)
dev.off()


pdf(paste0(currentdate,".NGSvsNoNGS_discoveries_per_year-allyears.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(3,4,1,0))
bars.x <- barplot(current.NGSvsnoNGS.allyears, beside=FALSE, ylab="Approximate # of gene discoveries by method", ylim=c(0,300), col=c(reviewred, reviewblue), pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=colnames(current.NGSvsnoNGS.allyears), cex=0.3, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.8, line=2)
# text(x=bars.x, y=colSums(current.NGSvsnoNGS.allyears), labels=colSums(current.NGSvsnoNGS.allyears), pos=3, xpd=TRUE, cex=0.3)
# text(x=bars.x[current.NGSvsnoNGS.allyears[2,]!=0], y=current.NGSvsnoNGS.allyears[1,]+0.5*current.NGSvsnoNGS.allyears[2,], labels=current.NGSvsnoNGS.allyears[2,current.NGSvsnoNGS.allyears[2,]!=0], xpd=TRUE, cex=0.5, col="white")
legend("topleft", fill=c(reviewred, reviewblue), legend=c("conventional", "NGS"), border=NA, bty="n", cex=0.8)
dev.off()

pdf(paste0(currentdate,".discoveries_per_year-cumulative-1986on.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(2.8,4,1,0))
plot(1:length(current.NGSvsnoNGS.since1986.cumsum), current.NGSvsnoNGS.since1986.cumsum, ylim=c(0,5000), type="o", pch=19, cex=0.65, cex.axis=0.7, cex.lab=0.8, xaxt="n", xlab="", ylab="Cumulative gene discoveries reported", bty="n")
mtext("Year", side=1, cex=0.7, line=1.9)
text(x=1:length(current.NGSvsnoNGS.since1986.cumsum), y=-700, labels=years.important, cex=0.7, xpd=TRUE, srt=90)
# axis(side=1, tick=TRUE, labels=FALSE, line=0.5)
points(1:length(current.NGSvsnoNGS.since1986.cumsum), cumsum(current.NGSvsnoNGS.since1986[1,]), type="b", col=reviewred, pch=19, cex=0.65)
lines(1:length(current.NGSvsnoNGS.since1986.cumsum), cumsum(current.NGSvsnoNGS.since1986[1,]), col=reviewred, pch=19)
points(1:length(current.NGSvsnoNGS.since1986.cumsum), cumsum(current.NGSvsnoNGS.since1986[2,]), type="b", col=reviewblue, pch=19, cex=0.65)
lines(1:length(current.NGSvsnoNGS.since1986.cumsum), cumsum(current.NGSvsnoNGS.since1986[2,]), col=reviewblue, pch=19)
legend("topleft", fill=c("black", reviewred, reviewblue), legend=c("total discoveries", "conventional", "NGS"), border=NA, bty="n", cex=0.8)
dev.off()




# graph rolling mean of discovery and syndrome delineation rates


baseR.rollmean <- function(dat, window) {
  n <- length(dat)
  y <- dat[window:n] - dat[c(1, 1:(n-window))]
  y[1] <- sum(dat[1:window])
  return(cumsum(y) / window)
}


current.NGSvsnoNGS.allyears.rollmean <- baseR.rollmean(colSums(current.NGSvsnoNGS.allyears), 3)
geneknownvsdelineation.allyears.rollmean <- baseR.rollmean(colSums(geneknownvsdelineation.allyears), 3)


pdf(paste0(currentdate,".discoveries_delineations_rollmean.pdf"), width=6, height=3.3, pointsize=12)
# par(mar=c(3,4,1,1))
par(mar=c(3,4,1,1))
plot(0, 0, xlim=c(1900,maxyear), ylim=c(0,300), xlab="", ylab="Annual total (rolling mean)", cex.axis=0.80, bty="n", xaxt="n")
	lines(names(current.NGSvsnoNGS.allyears.rollmean), current.NGSvsnoNGS.allyears.rollmean, col=reviewred, pch=21, cex=0.8)
	lines(names(geneknownvsdelineation.allyears.rollmean), geneknownvsdelineation.allyears.rollmean, col=phenoteal, pch=21, cex=0.8)
	lines(c("2012", "2015", "2018"), c(0, 120, 244), type="b", lty="dashed", pch=21, cex=0.8, col=cmgcolor)
	lines(1900:2011,rep(0, times=length(1900:2011)), lty="dashed", cex=0.8, col=cmgcolor)
	mtext("Year", side=1, cex=0.8, line=2)
	axis(side=1, at=seq(1900,2020,by=10), cex=0.8, tick=TRUE, line=0.5, labels=rep("",length(seq(1900,2020,by=10))))
	# axis(side=1, at=seq(1900,2020,by=10), cex=0.2, cex.lab=0.2, tick=TRUE, line=0.5, labels=seq(1900,2020,10))
	text(x=seq(1900,2020,by=10), y=-50, labels=seq(1900,2020,by=10), cex=0.8, xpd=TRUE)
	abline(v=maxyear, lty="dotted", col="#888888")
	text(x=maxyear, y=10, label=maxyear, pos=4, cex=0.7, col="#888888", xpd=TRUE)
	legend("topleft", fill=c(reviewred, phenoteal, cmgcolor, morecolors[7]), legend=c("reports of gene discoveries for MCs", "reports of new MCs", "gene discoveries for MCs by CMGs"), border=NA, bty="n", cex=0.8)
dev.off()

pdf(paste0(currentdate,".discoveries_delineations_rollmean_withDN.pdf"), width=6, height=3.3, pointsize=12)
par(mar=c(3,4,1,1))
plot(0, 0, xlim=c(1900,maxyear), ylim=c(0,300), xlab="", ylab="Annual total (rolling mean)", cex.axis=0.80, bty="n", xaxt="n")
	lines(names(current.NGSvsnoNGS.allyears.rollmean), current.NGSvsnoNGS.allyears.rollmean, col=reviewred, pch=21, cex=0.8)
	lines(names(geneknownvsdelineation.allyears.rollmean), geneknownvsdelineation.allyears.rollmean, col=phenoteal, pch=21, cex=0.8)
	lines(names(solved.allmodels["de novo",]), solved.allmodels["de novo",], col=morecolors[7], cex=0.8, lty="dotted")
	lines(c("2012", "2015", "2018"), c(0, 120, 244), type="b", lty="dashed", pch=21, cex=0.8, col=cmgcolor)
	lines(1900:2011,rep(0, times=length(1900:2011)), lty="dashed", cex=0.8, col=cmgcolor)
	mtext("Year", side=1, cex=0.8, line=2)
	axis(side=1, at=seq(1900,2020,by=10), cex=0.8, tick=TRUE, line=0.5, labels=rep("",length(seq(1900,2020,by=10))))
	text(x=seq(1900,2020,by=10), y=-50, labels=seq(1900,2020,by=10), cex=0.8, xpd=TRUE)
	abline(v=maxyear, lty="dotted", col="#888888")
	text(x=maxyear, y=10, label=maxyear, pos=4, cex=0.7, col="#888888", xpd=TRUE)
	legend("topleft", fill=c(reviewred, phenoteal, cmgcolor, morecolors[7]), legend=c("reports of gene discoveries for MCs", "reports of new MCs", "gene discoveries for MCs by CMGs", "reports of MCs caused by DNVs"), border=NA, bty="n", cex=0.8)
dev.off()











pdf(paste0(currentdate,".discoveries_per_year_inheritance.pdf"), width=6, height=4, pointsize=12)
par(mar=c(3,4,3,0))
bars.x <- barplot(solved.allmodels, beside=FALSE, ylab="Approximate # of gene discoveries", ylim=c(0,300), col=mycolors, pch=21, border=NA, cex.main=0.95, cex.axis=0.8, xaxt="n")
text(x=bars.x, y=-20, labels=years.important, cex=0.8, xpd=TRUE, srt=90)
mtext("Year", side=1, cex=0.8, line=2)
# text(x=bars.x, y=colSums(solved.allmodels), labels=colSums(solved.allmodels), pos=3, xpd=TRUE, cex=0.5)
# text(x=bars.x[solved.allmodels[2,]!=0], y=solved.allmodels[1,]+0.5*solved.allmodels[2,], labels=current.NGSvsnoNGS.since1986[2,current.NGSvsnoNGS.since1986[2,]!=0], xpd=TRUE, cex=0.5, col="white")
legend("topleft", fill=rev(mycolors), legend=rev(modelnames), border=NA, bty="n", cex=0.8)
dev.off()
