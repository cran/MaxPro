\name{MaxProAugment}
\alias{MaxProAugment}

\title{
Augment an Existing Design based on the Maximum Projection Criterion
}
\description{
Select the best set of design points to augment a given design matrix by optimizing the maximum projection (MaxPro) criterion sequentially. When multiple types of factors exist, all continuous factors should be placed as the leftmost columns in the design matrix, discrete numeric factors should be placed as middle columns, and all nominal factors should be placed at the rightmost columns in the design matrix. The ordinal factors, if any, should be pre-converted into discrete numeric factors through the scoring method (see, e.g., Wu and Hamada 2009, Section 14.10). This function can also be used for creating a nonadaptive sequential design, designs to run in batches, validation runs, designs in constrained regions, and nested designs.
}
\usage{
MaxProAugment(ExistDesign, CandDesign, nNew, 
p_disnum=0, l_disnum=NULL, p_nom=0, l_nom=NULL)
}

\arguments{
  \item{ExistDesign}{
The existing design matrix to be augmented, where each row is an experimental run and each column is a factor. All columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1].
}
  \item{CandDesign}{
The matrix containing candidate design points, whose columns should match those in the ExistDesign matrix. All columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1]. Function \code{\link{CandPoints}} can be used to quickly generate candidate design points.
}
  \item{nNew}{
The number of new design points to add 
}
  \item{p_disnum}{
Optional, default is 0. The number of discrete numeric factors
}
  \item{l_disnum}{
Optional, default is a vector containing the number of distinct levels for each discrete numeric factor in the ExistDesign and CandDesign. The vector specifying the number of levels for each discrete numeric factor 
}
  \item{p_nom}{
Optional, default is 0. The number of nominal factors
}
  \item{l_nom}{
Optional, default is a vector containing the number of levels for each nominal factor in the ExistDesign and CandDesign. The vector specifying the number of levels for each nominal factor 
}
}
\details{
This function uses a greedy search algorithm to select the best nNew points out of the CandDesign matrix to augment the ExistDesign matrix based on the MaxPro criterion. Factor constraints (irregular-shaped design space), if any, can be incorporated through the specification of valid candidate design points in CandDesign (please see the example below). 
}

\value{
The value returned from the function is a list containing the following components:
\item{Design}{The augmented design matrix (including the original existing design points)}
\item{measure}{The MaxPro criterion measure of the augmented design}
\item{time_rec}{Time to complete the search}
}
\references{
Joseph, V. R., Gul, E., and Ba, S. (2015) "Maximum Projection Designs for Computer Experiments," \emph{Biometrika}, 102, 371-380.

Joseph, V. R. (2016) "Rejoinder," \emph{Quality Engineering}, 28, 42-44.

Joseph, V. R., Gul, E., and Ba, S. (2018) "Designing Computer Experiments with Multiple Types of Factors: The MaxPro Approach," \emph{Journal of Quality Technology}, to appear.

Wu, C. F. J., and Hamada, M. (2009), \emph{Experiments: Planning, Analysis, and Parameter Design Optimization, 2nd Edition}, New York: Wiley.
}
\author{
Shan Ba <shanbatr@gmail.com> and V. Roshan Joseph <roshan@isye.gatech.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CandPoints}}, \code{\link{MaxProLHD}}, \code{\link{MaxProQQ}}
}

\examples{

ExistDesign=matrix(c(0.11,0.96,0.42,0.99),ncol=2,byrow=TRUE)

#Generate candidate design points that satisfy the factor constraints
S=CandPoints(N=10000,p_cont=2)
S=S[S[,1]+S[,2]<=1.5,] #constraint 1
S=S[(S[,1]-.5)^2+(S[,2]-.5)^2>=(2*.25^2),] #constraint 2
plot(S[,1],S[,2],bty="n",cex=.5,xlab=expression(x[1]),ylab=expression(x[2]),xlim=c(0,1),ylim=c(0,1))

result=MaxProAugment(ExistDesign, CandDesign=S, nNew=20)
D=result$Design
points(D[,1],D[,2],col='red',pch=16)
text(D[,1],D[,2],1:nrow(D),col=2,cex=2)


}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
