\name{MaxProImproveLHD}
\alias{MaxProImproveLHD}

\title{
Improve a Latin Hypercube Design based on the Maximum Projection Criterion
}
\description{
Improve an initial Latin hypercube design by permuting its levels for each variable based on the maximum projection (MaxPro) criteion
}
\usage{
MaxProImproveLHD(InitialDesign,s=2,localopm=FALSE,temp0=0,nstarts=1,itermax=400,
total_iter=1e+06)
}

\arguments{
  \item{InitialDesign}{
The initial Latin Hypercube design matrix to be improved
}
  \item{s}{
Optional, default is ``2''. The parameter in defining the s-norm distance (2 corresponds to Euclidean distance)
}
  \item{localopm}{
Optional, default is "FALSE" and the algorithm conducts global searches using simulated annealing; if "TURE", the algorithm only conducts local searches
}
  \item{temp0}{
Optional, The initial temperature in the simulated annealing algorithm. Change this value if you want to start with a higher or lower temperature 
}
  \item{nstarts}{
Optional, default is ``1''. The number of random starts
}
  \item{itermax}{
Optional, default is ``400''. The maximum number of non-improving searches allowed under each temperature. Lower this parameter if you want the algorithm to converge faster 
}
  \item{total_iter}{
Optional, default is ``1e+06''.The maximum total number of iterations. Lower this number if the design is prohibitively large and you want to terminate the algorithm prematurely to report the best design found so far
}
}
\details{
This function utilizes a version of the simulated annealing algorithm to improve a given Latin hypercube design by permuting its levels for each variable based on the MaxPro criterion. Parameters in the algorithm may need to be properly tuned to achieve global convergence. Please refer to Joseph, Gul and Ba (2015) for details.
}

\value{
The value returned from the function is a list containing the following components:
\item{Design}{The improved Latin hypercube design matrix}
\item{temp0}{Initial temperature for simulated annealing}
\item{measure}{The MaxPro criterion measure}
\item{time_rec}{Time to complete the search}
\item{ntotal}{The total number of iterations}
}
\references{
Joseph, V. R., Gul, E., and Ba, S. (2015) "Maximum Projection Designs for Computer Experiments," \emph{Biometrika}.
}
\author{
Shan Ba <shanbatr@gmail.com> and V. Roshan Joseph <roshan@isye.gatech.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

InitialDesign<-matrix(rep(seq(0,1,length=10),4),ncol=4)
DOX<-MaxProImproveLHD(InitialDesign)
DOX$Design

}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
