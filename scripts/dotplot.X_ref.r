args = commandArgs(T)

#========== for the dot plot ==========

#plotFilenameTemplate = "hg38Y_onto_hg002Y.dotplot.pdf"
plotFilenameTemplate = args[5]

addTitle = F
addLegend = F
delineateChromosome = F

#~~~~~ read the hg002 and hg38 lengths

#scaffolds.hg002 = read_scaffold_lengths("chrY_hg002_v2.7.lengths")
#xLength = sum(scaffolds.hg002$length)
xLength = read.table(args[1], h=F)
xLength = xLength[1,2]

#scaffolds.hg38 = read_scaffold_lengths("hg38.chrY.lengths",sortByLength=F)
#yLength = sum(scaffolds.hg38$length)
yLength = read.table(args[2], h=F)
yLength = yLength[1,2]
#~~~~~ read the dotplot

colClasses=c("character","integer","character","integer","integer")
#dots = read.table("hg38Y_onto_hg002Y.dotplot.gz",header=T,comment.ch="",colClasses=colClasses)
dots = read.table(args[3], header=T,comment.ch="")

#~~~~~ set up hg002 region boundaries
# these names and colors manually taken from Monika_HG002_classes.v5.dat

regionNames = c("PAR","X-DEG","XTR","AMPL","SAT","CEN","DYZ","HET","OTHER")

regionColor = c(rgb(151,203,153,max=255), # PAR
				rgb(255,239,87,max=255),  # X-DEG
				rgb(238,169,186,max=255), # XTR
				rgb(136,192,234,max=255), # AMPL
				rgb(106,71,0,max=255),    # SAT
				rgb(176,32,38,max=255),   # CEN
				rgb(106,71,0,max=255),    # DYZ
				rgb(119,119,119,max=255), # HET
				rgb(217,216,216,max=255)) # OTHER
names(regionColor) = regionNames

#colClasses=c("character","numeric","numeric","character","numeric","character","numeric","numeric","character")
#regions = read.table("Monika_HG002_classes.v5.dat",header=T,comment.ch="",colClasses=colClasses)
colClasses=c("character","numeric","numeric","character")
regions = read.table(args[4], header=T,comment.ch="",colClasses=colClasses)

regions$class = gsub("^PAR.*$","PAR",regions$class)
regions$class = gsub("^XTR.*$","XTR",regions$class)
regions$class = gsub("^DYZ.*$","DYZ",regions$class)

legText  = regionNames
legColor = regionColor[legText]

#~~~~~ plot

xSpan = xLength
ySpan = yLength
xlim  = c(0,xLength)
ylim  = c(0,yLength)

windowSize = 12
aspect = 1.069 # this makes the plot area ≈ square if genomes are same length
aspect = aspect*ySpan/xSpan
width  = ifelse(aspect<=1,windowSize,windowSize/aspect);
height = width*aspect;
height = height*0.9909365    # special adjustment to match horizontal axis

plotFilename = NULL
if (!is.null(plotFilenameTemplate))
	plotFilename = plotFilenameTemplate

turnDeviceOff = F

if (is.null(plotFilename)) {
	quartz(width=width,height=height)
} else {
	print(paste("drawing to",plotFilename))
	pdf(file=plotFilename,width=width,height=height,pointsize=9)
	turnDeviceOff = T
}

#title = "lastz alignment, HG002 Y vs MF2 Y"
#xlab  = "HG002 Y"
#ylab  = "MF2 Y"
xlab = basename(args[1])
ylab = basename(args[2])

# create window
options(scipen=10)
par(mar=c(4,4.1,2.5,0.2)+0.1);    # BLTR
if (addTitle)
{
	plot(NA,xlim=xlim,ylim=ylim,main=title,xlab=xlab,ylab=ylab)
} else {
	plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
}

# draw background regions
for (ix in 1:(nrow(regions)))
{
	xLeft   = regions$start[ix]
	xRight  = regions$end[ix]
	color   = regionColor[regions$class[ix]]
	rect(xLeft,0,xRight,yLength,border=NA,col=color)
}

# plot the dotplot
#lines(dots$pos1,dots$pos2,lwd=1,col="black")
#segments(dots$zstart2., dots$zstart1, dots$end2., dots$end1, lwd=1,col="black")
segments(dots$zstart1, dots$zstart2., dots$end1, dots$end2., lwd=1,col="black")

# delineate the chromosome
if (delineateChromosome)
{
	lines(c(0,xLength),c(0,0),col="red")
	lines(c(0,xLength),c(yLength,yLength),col="red")
	lines(c(0,0),c(0,yLength),col="red")
	lines(c(xLength,xLength),c(0,yLength),col="red")
}

if (addLegend)
	legend("bottomright",bg="white",cex=1.0,legText,fill=legColor)

if (turnDeviceOff) dev.off()

