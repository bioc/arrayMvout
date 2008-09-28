\name{arrOutStruct-class}
\docType{class}
\alias{arrOutStruct-class}
\alias{plot,arrOutStruct,ANY-method}
\alias{show,arrOutStruct-method}

\title{Class "arrOutStruct" container for ArrayOutliers output }
\description{Class "arrOutStruct" container for ArrayOutliers output }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("arrOutStruct", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} includes
components outl, inl, with results of outlier detection applied
to principal components of multivariate QA statistics }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "arrOutStruct", y = "ANY")}: a biplot of QA statistics }
    \item{show}{\code{signature(object = "arrOutStruct")}: summary report }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Vince Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("arrOutStruct")
}
\keyword{classes}
