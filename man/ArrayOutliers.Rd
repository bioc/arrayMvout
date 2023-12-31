\name{ArrayOutliers}
\alias{ArrayOutliers}
\alias{afxsubDEG}
\alias{afxsubQC} 
\alias{s12c} 
\alias{s12cDEG} 
\alias{s12cQC}
\alias{ILM1}
\alias{itnQA}
\alias{maqcQA}
\alias{ilmQA}
\alias{spikQA}
\alias{fig3map}
\title{Multivariate outlier detection based on PCA of
QA statistics  }
\description{Multivariate outlier detection based on PCA of
QA statistics  }
\usage{
ArrayOutliers (data, alpha, alphaSeq = c(0.01, 0.05, 0.1),  ... )
#    qcOutput = NULL, plmOutput = NULL, degOutput = NULL, prscale = TRUE, 
#    pc2use = 1:3) 
}
\arguments{
  \item{data}{ an (affy) AffyBatch instance with at least 11 samples}
  \item{alpha}{ false positive rate for outlier detection, adjusting
for multiple comparisons according to Caroni and Prescott's
adaptation of Rosner (1983); full report based on this choice of alpha}
  \item{alphaSeq}{vector of alpha candidates to be quickly tried for
short report}
  \item{\dots}{ additional parameters, see below }
}
\details{
Additional parameters may be supplied
\describe{
  \item{qcOutput}{optional result of simpleaffy qc() to speed computations}
  \item{plmOutput}{optional result of affyPLM fitPLM() to speed computations}
  \item{degOutput}{optional result of affy AffyRNAdeg() to speed computations}
  \item{prscale}{scaling option for prcomp}
  \item{pc2use}{selection of principal components to use for outlier detection}
}
Data elements afxsubDEG, afxsubQC, s12cDEG, s12cQC are precomputed RNA degradation
and simpleaffy qc() results; s12c is an AffyBatch with digital contamination of some
samples.

Data elements maqcQA and itnQA are affymetrix QC statistics on large
collections of arrays.  Data element ilmQA is a derived from
a LumiBatch of the Illumina-submitted
MAQC raw data, 19 arrays.  (Conveyed by Leming Shi, personal communication).
Data element spikQA is a 12x9 matrix of QA parameters obtained for
12 arrays from U133A spikein dataset, with first 2 arrays digitally contaminated
as described in  Asare et al.

Data element fig3map gives the indices of the points labeled A-H in
Figure 3 of the manuscript by Asare et al. associated with this package.
}
\value{
an instance of arrOutStruct class, a list with a partition of samples into two data frames (inl and
outl) with QA summary statistics 
}
%\references{  }
\author{Z. Gao et al.}
%\note{  }


%\seealso{  }

\examples{
\dontrun{
library(simpleaffy)
setQCEnvironment("hgu133acdf")  # no CDF corresponding to tag array
if ( require("mvoutData") ) { 
 data(s12c)
 data(s12cQC)
 data(s12cDEG)
 library(affyPLM)
 s12cPset = fitPLM(s12c)
 ao = ArrayOutliers(s12c, alpha=0.05, qcOut=s12cQC, plmOut=s12cPset, degOut=s12cDEG)
 ao
 }
if (require("lumiBarnes")) {
 library(lumiBarnes)
 data(lumiBarnes)
 ArrayOutliers(lumiBarnes, alpha=0.05)
 lb2 = lumiBarnes
 exprs(lb2)[1:20000,1:2] = 10000*exprs(lb2)[1:20000,1:2]
 ArrayOutliers(lb2, alpha=0.05)
 }
data(maqcQA) # affy
ArrayOutliers(maqcQA[,-c(1:2)], alpha=.05)
ArrayOutliers(maqcQA[,-c(1:2)], alpha=.01)
data(ilmQA) # illumina
ArrayOutliers(data.frame(ilmQA), alpha=.01)
data(itnQA) # 507 arrays from ITN
ArrayOutliers(itnQA, alpha=.01)
}
}
\keyword{ models }
