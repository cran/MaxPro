\name{MaxProRunOrder}
\alias{MaxProRunOrder}

\title{
Find the Optimal Sequential Order to Run a Given Experimental Design
}
\description{
Given a computer experimental design matrix, this function searches for an optimal run (row) order based on the maximum projection (MaxPro) criterion. This optimal order enables the given design to be run in a sequential manner: when terminated at any step, the previous design points form a nearly optimal subset based on the MaxPro criterion. 
}
\usage{
MaxProRunOrder(Design,p_nom=0,initial_row=1)
}

\arguments{
  \item{Design}{
The design matrix, where each row is an experimental run and each column is a factor. The rightmost p_nom columns correspond to the p_nom nominal factors, and the columns on the left are for continuous factors and discrete numeric factors. The ordinal factors, if any, should be pre-converted into discrete numeric factors through the scoring method (see, e.g., Wu and Hamada 2009, Section 14.10). All columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1].
}
  \item{p_nom}{
Optional, default is 0. The number of nominal factors
}
  \item{initial_row}{
Optional, default is 1. The vector specifying the row number of each design point in the given design matrix that should be run at first or have already been run. 
}
}
\details{
This function utilizes a greedy search algorithm to find the optimal row order to run the given experimental design based on the MaxPro criterion.  
}

\value{
The value returned from the function is a list containing the following components:
\item{Design}{The design matrix in optimal run (row) order. The run sequence ID is added as the first column}
\item{measure}{The MaxPro criterion measure of the given design}
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
\code{\link{MaxProLHD}}, \code{\link{MaxProQQ}}, \code{\link{MaxProAugment}}
}

\examples{
D0=MaxProLHD(25,2)$Design 

#Assume the first two rows of the design have already been executed 
#Find the optimal run orders

D=MaxProRunOrder(D0,p_nom=0,initial_row=c(1,2))$Design
plot(D[,2],D[,3],xlim=c(0,1),ylim=c(0,1),type="n",
xlab=expression(x[1]),ylab=expression(x[2]),cex.lab=1.5)
text(D[,2],D[,3],labels=D[,1],col='red')
}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
