\name{ArrayOutliers-methods}
\docType{methods}
\alias{ArrayOutliers-methods}
\alias{ArrayOutliers,ANY,missing,missing-method}
\alias{ArrayOutliers,AffyBatch,numeric,ANY-method}
\alias{ArrayOutliers,LumiBatch,numeric,ANY-method}
\title{ ArrayOutliers -- wrapper for platform-specific multivariate outlier detection 
 for expression arrays}
\description{
wraps functions that perform multivariate outlier detection on dimension-reduced QA statistics
of expression arrays
}
\section{Methods}{
\describe{

\item{data = "ANY", alpha = "missing", alphaSeq = "missing"}{ fails; tells user that alpha is obligatory
parameter }

\item{data = "AffyBatch", alpha = "numeric", alphaSeq = "ANY"}{ performs calibrated
multivariate outlier detection on an AffyBatch instance using various affy-specific QA
parameters }

\item{data = "LumiBatch", alpha = "numeric", alphaSeq = "ANY"}{ performs calibrated
multivariate outlier detection on an LumiBatch instance using various illumina-specific QA
parameters }
}}
\examples{
example(ArrayOutliers)
}
\keyword{methods}
