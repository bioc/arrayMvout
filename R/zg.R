setClass("arrOutStruct", representation(call="call"), contains="list")

###################################
### Microarray QA outlier detection
### Jan 8, 2008
###################################

#########################################################
### Function for microarray outlier detection
### Input: AffyBatch object
### Output: A list of QA outlier arrays with QA metrics
###         Non-outlier arrays with QA metrics
#########################################################

setGeneric("ArrayOutliers", function(data, alpha, alphaSeq=c(0.01, 0.05, 0.10), ...) {
 standardGeneric("ArrayOutliers") })

setMethod("ArrayOutliers", c("ANY", "missing", "missing"), 
   function( data, alpha, alphaSeq ) stop("must supply second argument: false outlier labeling rate"))
   
 
setMethod("ArrayOutliers", c("AffyBatch", "numeric", "ANY"),
  function(data, alpha, alphaSeq=c(.01, .05, .10), 
   qcOutput=NULL, plmOutput=NULL, degOutput=NULL, prscale=TRUE, 
   pc2use = 1:3){
   .affyArrayOutliers( data, alpha, alphaSeq, qcOutput, plmOutput,
     degOutput, prscale, pc2use ) } )

.affyArrayOutliers <- function(data, alpha=.05, alphaSeq=c(0.01, 0.05, 0.10),
   qcOutput=NULL, plmOutput=NULL, degOutput=NULL, prscale=TRUE, 
   pc2use = 1:3){

#################################
### QA metrics
#################################

# VJC -- I have lots of problems with names and list operations
# used in the original source -- 

pd = data.frame(pData(data));
sn = sampleNames(data)

### Summary of Affymetrix QC metrics
if (is.null(qcOutput)) Data.qc <- qc(data)
else Data.qc = qcOutput

Affy.QC.Summary <- data.frame(cbind(avbg(Data.qc),
				    sfs(Data.qc ),
				    percent.present(Data.qc),
				    ratios(Data.qc)[, 1:2]));
names(Affy.QC.Summary)<- c("avgBG", "SF", "Present", "HSACO7", "GAPDH");

Affy.QC.Summary$samp <- sn

### affyPLM
if (is.null(plmOutput)) Pset <- fitPLM(data)
else Pset = plmOutput

nuse <- boxplot(Pset, plot=FALSE);

rle <- Mbox(Pset, plot=FALSE);

z <- t(rbind(nuse$stats, rle$stats ) );
PLM <- data.frame(z);
names(PLM) <- c( "nuse.lw", "nuse.q25", "NUSE", "nuse.q75", "nuse.uw", "rle.lw", "rle.q25", "RLE", "rle.q75", "rle.uw" );
PLM$RLE_IQR <- PLM$rle.q75-PLM$rle.q25;
PLM$samp = sn

### Compute RNA degrdation slope
if (is.null(degOutput)) RNAdeg <- AffyRNAdeg(data)
else RNAdeg = degOutput

RNAdegSummary <- data.frame(samp=sn, slope=RNAdeg$slope, 
   pvalue=RNAdeg$pvalue);

# check for compatibility
if (!(all.equal(Affy.QC.Summary$samp, PLM$samp))) stop("affyqc df and plm df out of synch")
if (!(all.equal(PLM$samp, as.character(RNAdegSummary$samp)))) stop("rnadeg df and plm df out of synch")

### merge Affy, PLM, and RNAdeg slop -- NO!  avoid merge
#AffyPLM <- merge(Affy.QC.Summary, PLM, by="samp");
#AffyPLMRNA <- merge(AffyPLM, RNAdegSummary, by="samp");

######################################
### QA Outlier detection
######################################
QA <- data.frame(samp=PLM$samp,  # confirmed to be in right order
                 avgBG=Affy.QC.Summary$avgBG,
                 SF=Affy.QC.Summary$SF,
                 Present=Affy.QC.Summary$Present,
                 HSACO7=Affy.QC.Summary$HSACO7,
                 GAPDH=Affy.QC.Summary$GAPDH,
                 NUSE=PLM$NUSE,
                 RLE=PLM$RLE,
                 RLE_IQR=PLM$RLE_IQR,
                 RNAslope=RNAdegSummary$slope)

### back transform log2-based HSACO7 and GAPDH
QAback = QA
QAback[, "HSACO7"] = 2^QAback[, "HSACO7"]
QAback[, "GAPDH"] = 2^QAback[, "GAPDH"]


### PCA
princomp.matrix <- QA[ , c(2:10)];
#princomp.object <- princomp(princomp.matrix, scale=prscale, cor=TRUE, scores=TRUE);
princomp.object <- prcomp(princomp.matrix, scale=prscale)

### PCA followed by Sequential Wilk's test 

QCPC = princomp.object$x[, pc2use]

PC.wilk.all <- lapply(alphaSeq, function(alpha)mv.calout.detect(QCPC, alpha=alpha));
PC.wilk.results <- mv.calout.detect(QCPC, alpha=alpha)

### Outlier list
QCPCwilk <- PC.wilk.results$inds;
QCPCoutlier <- NULL
if (length(QCPCwilk)>0) QCPCoutlier = sn[QCPCwilk]

### data with outliers excluded
QCPCNoOutliers <- subset(QCPC, !is.element(sn, QCPCoutlier));


### QA profiles of outliers
QCoutlier <- subset(QAback, is.element(QAback$samp, QCPCoutlier));
#write.table(QCoutlier, 'Array_QA_outliers.csv', col.names=NA, sep=",");

### QA profile with outliers excluded
QCnooutlier <- subset(QAback, !is.element(QAback$samp, QCPCoutlier));
#write.table(QCnooutlier, 'Array_QA_outlier_Excluded.csv', col.names=NA, sep=",");
inl = QCnooutlier
rownames(inl) = inl$samp
inl = inl[-1]
rownames(QAback)  = QAback$samp
QAback = QAback[-1]
ans = list(outl=QCoutlier, inl=inl, QA=QAback, PC.wilk.all=PC.wilk.all, 
 alphaSeq=alphaSeq,
 Pset=Pset, QCset=Data.qc, DEGset=RNAdeg)

new("arrOutStruct", call=match.call(), ans)
}

setMethod("show", "arrOutStruct", function(object) {
cat("ArrayOutliers result.\n")
if (is.null(object[[1]])) {
 cat("No outliers. First row of QC features\n")
 print(object[[3]][1,,drop=FALSE])
}
else {
  cat("There were ", nsamp <- nrow(object[[3]]), " samples with ",
   nout <- nrow(object[[1]]), " outliers detected.\n")
  if (nout > 0) {
   cat("Coordinate-wise means of inlying arrays:\n")
   print(apply(object$inl,2,mean,na.rm=TRUE))
   cat("Features of outlying arrays:\n")
   print(object$outl)
 }
}
})
  
   
setMethod("plot", "arrOutStruct", function(x, y, ..., main="all QC stats") {
 biplot(prcomp(x$QA),main=main, ...)
})


setMethod("ArrayOutliers", c("LumiBatch", "numeric", "ANY"),
  function(data, alpha, alphaSeq=c(.01, .05, .10), 
   qcOutput=NULL, plmOutput=NULL, degOutput=NULL, prscale=TRUE, 
   pc2use = 1:3){
   .lumiArrayOutliers( data, alpha, alphaSeq, qcOutput, plmOutput,
     degOutput, prscale, pc2use ) } )

.lumiArrayOutliers <- function(data, alpha=.05, alphaSeq=c(0.01, 0.05, 0.10),
   qcOutput=NULL, plmOutput=NULL, degOutput=NULL, prscale=TRUE, 
   pc2use = 1:3){
   LQ = t(lumiQ(data)@QC$sampleSummary)
   PLQ = prcomp(LQ, scale=prscale)$x[, pc2use]
   odat = mv.calout.detect(PLQ,alpha=alpha)
   PC.wilk.all <- lapply(alphaSeq, function(alpha)mv.calout.detect(PLQ, alpha=alpha));
   if (is.na(odat$inds[1])) {
       QCoutlier=NULL
       QCnooutlier=LQ
       }
   else {
       QCoutlier=LQ[odat$inds,,drop=FALSE]
       QCnooutlier=LQ[-odat$inds,,drop=FALSE]
       }
ans = list(outl=QCoutlier, inl=QCnooutlier, QA=PLQ, PC.wilk.all=PC.wilk.all, 
 alphaSeq=alphaSeq, Pset=NULL, QCset=LQ, DEGset=NULL)
new("arrOutStruct", call=match.call(), ans)
}
   
   
setMethod("ArrayOutliers", c("data.frame", "numeric", "ANY"),
  function(data, alpha, alphaSeq=c(.01, .05, .10), pc2use=1:3) {
  pp = prcomp(data, scale=TRUE)$x[,pc2use]
  odat = mv.calout.detect(pp, alpha=alpha)
  pall = lapply(alphaSeq, function(alpha)mv.calout.detect(pp, alpha=alpha))
   if (is.na(odat$inds[1])) {
       QCoutlier=NULL
       QCnooutlier=data
       }
   else {
       QCoutlier=data[odat$inds,,drop=FALSE]
       QCnooutlier=data[-odat$inds,,drop=FALSE]
       }
ans = list(outl=QCoutlier, inl=QCnooutlier, QA=data, PC.wilk.all=pall,
 alphaSeq=alphaSeq, Pset=NULL, QCset=NULL, DEGset=NULL)
new("arrOutStruct", call=match.call(), ans)
})
 
