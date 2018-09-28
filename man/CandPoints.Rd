\name{CandPoints}
\alias{CandPoints}

\title{
Generate Candidate Design Points Randomly for Various Types of Factors
}
\description{
Quickly generate random candidate design points for different types of factors which can be used as an input to the \code{\link{MaxProAugment}} function.
}
\usage{
CandPoints(N,p_cont,l_disnum=NULL,l_nom=NULL)
}

\arguments{
  \item{N}{
The number of random candidate design points to generate
}
  \item{p_cont}{
The number of continuous factors
}
  \item{l_disnum}{
Optional, default is ``NULL''. A vector containing the number of levels for each discrete numeric factor (NULL if there are no discrete numeric factors)
}
  \item{l_nom}{
Optional, default is ``NULL''. A vector containing the number of levels for each nominal factor (NULL if there are no nominal factors)
}
}

\details{
This function generates uniform random numbers for each continuous factor and randomly sampled levels for each discrete numeric and nominal factor. The result can be used as the candidate design in the \code{\link{MaxProAugment}} function.
}

\value{
The value returned from the function is a matrix containing the candidate design points. 
}

\author{
Shan Ba <shanbatr@gmail.com> and V. Roshan Joseph <roshan@isye.gatech.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MaxProAugment}}
}

\examples{

CCC<-CandPoints(N=10000,p_cont=3,l_disnum=c(3,3,4),l_nom=c(2,2,3))

}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
