\name{MaxPro-package}
\alias{MaxPro-package}
\docType{package}
\title{
MaxPro package
}
\description{
Maximum projection (MaxPro) designs for quantitative and/or qualitative factors
}
\details{
\tabular{ll}{
Package: \tab MaxPro\cr
Type: \tab Package\cr
Version: \tab 4.1-2\cr
Date: \tab 2018-09-27\cr
License: \tab LGPL-2.1\cr
}
Important functions in this package are: \code{\link{MaxProLHD}} generates the MaxPro Latin hypercube design for continuous factors, \code{\link{MaxPro}} searches for the locally optimal MaxPro design in the neighborhood of a given initial design for continuous factors, \code{\link{MaxProQQ}} generates the MaxPro design with both quantitative and qualitative factors, \code{\link{MaxProAugment}} sequentially augments an existing design with additional design points based on the MaxPro criterion for all types of factors and constraints, \code{\link{MaxProRunOrder}} finds the optimal order to run a given design in a sequential manner based on the MaxPro criterion, and \code{\link{MaxProMeasure}} reports the MaxPro criterion value for a given design.
}
\author{
Shan Ba and V. Roshan Joseph

Maintainer: Shan Ba <shanbatr@gmail.com>
}
\references{
Joseph, V. R., Gul, E., and Ba, S. (2015) "Maximum Projection Designs for Computer Experiments," \emph{Biometrika}, 102, 371-380.

Joseph, V. R., Gul, E., and Ba, S. (2018) "Designing Computer Experiments with Multiple Types of Factors: The MaxPro Approach," \emph{Journal of Quality Technology}, to appear.
}

\keyword{ package }

