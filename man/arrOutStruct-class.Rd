\name{arrOutStruct-class}
\docType{class}
\alias{arrOutStruct-class}
\alias{plot,arrOutStruct,ANY-method}
\alias{show,arrOutStruct-method}

\title{Class "arrOutStruct" container for ArrayOutliers output }
\description{Class "arrOutStruct" container for ArrayOutliers output }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("arrOutStruct", ...)}.
This class just extends \code{list} but has specialized \code{show}
and \code{plot} methods.
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{\link[Biobase:class.AssayData]{AssayData}}, by class "list", 
distance 2.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "arrOutStruct", y = "ANY")}: a biplot of QA statistics }
    \item{show}{\code{signature(object = "arrOutStruct")}: summary report }
	 }
}
\author{ Vince Carey <stvjc@channing.harvard.edu> }

\examples{
data(maqcQA)
f1 = ArrayOutliers(maqcQA[,-c(1:2)], alpha=0.01)
names(f1)
f1
}
\keyword{classes}
