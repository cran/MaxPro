\name{MaxProMeasure}
\alias{MaxProMeasure}

\title{
Maximum Projection Criterion Measure
}
\description{
Calculate the maximum projection (MaxPro) criterion measure for a given design. 
}
\usage{
MaxProMeasure(Design, p_nom=0)
}

\arguments{
  \item{Design}{
The design matrix where each row is an experimental run and each column is a factor. The rightmost p_nom columns correspond to the p_nom nominal factors, and the columns on the left are for continuous factors and discrete numeric factors. The ordinal factors, if any, should be pre-converted into discrete numeric factors through the scoring method (see, e.g., Wu and Hamada 2009, Section 14.10). All columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1].
}
  \item{p_nom}{
Optional, default is ``0''. The number of nominal factors
}
}
\details{
This function computes the MaxPro criterion measure for a given design. Please refer to Joseph, Gul and Ba (2015) and Joseph, Gul and Ba (2018) for details.
}

\value{
The value returned from the function is the MaxPro criterion measure
}
\references{
Joseph, V. R., Gul, E., and Ba, S. (2015) "Maximum Projection Designs for Computer Experiments," \emph{Biometrika}, 102, 371-380.

Joseph, V. R., Gul, E., and Ba, S. (2018) "Designing Computer Experiments with Multiple Types of Factors: The MaxPro Approach," \emph{Journal of Quality Technology}, to appear.

Wu, C. F. J., and Hamada, M. (2009), \emph{Experiments: Planning, Analysis, and Parameter Design Optimization, 2nd Edition}, New York: Wiley.
}
\author{
Shan Ba <shanbatr@gmail.com> and V. Roshan Joseph <roshan@isye.gatech.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MaxProLHD}}, \code{\link{MaxProQQ}}, \code{\link{MaxProAugment}}
}

\examples{

D<-MaxProLHD(10,2)$Design
MaxProMeasure(D,p_nom=0)


}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
