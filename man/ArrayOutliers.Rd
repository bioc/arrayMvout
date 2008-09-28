\name{ArrayOutliers}
\alias{ArrayOutliers}
\alias{afxsubDEG}
\alias{afxsubQC} 
\alias{s12c} 
\alias{s12cDEG} 
\alias{s12cQC}
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
library(simpleaffy)
setQCEnvironment("hgu133acdf")  # tag array not known...
data(s12c)
data(s12cQC)
data(s12cDEG)
library(affyPLM)
s12cPset = fitPLM(s12c)
ao = ArrayOutliers(s12c, alpha=0.05, qcOut=s12cQC, plmOut=s12cPset, degOut=s12cDEG)
ao
library(lumiBarnes)
data(lumiBarnes)
ArrayOutliers(lumiBarnes, alpha=0.05)
lb2 = lumiBarnes
exprs(lb2)[1:20000,1:2] = 10000*exprs(lb2)[1:20000,1:2]
ArrayOutliers(lb2, alpha=0.05)
}
\keyword{ models }
