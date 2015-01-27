\name{MaxProLHD}
\alias{MaxProLHD}

\title{
Maximum Projection Latin Hypercube Designs
}
\description{
Generate the maximum projection (MaxPro) Latin hypercube design based on a simulated annealing algorithm
}
\usage{
MaxProLHD(n, p, s=2, temp0=0, nstarts = 1, itermax = 400, total_iter = 1e+06)
}

\arguments{
  \item{n}{
The number of runs (design points)
}
  \item{p}{
The number of input factors (variables)
}
  \item{s}{
Optional, default is ``2''. The parameter in defining the s-norm distance (2 corresponds to Euclidean distance)
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
This function utilizes a version of the simulated annealing algorithm to efficiently generate the optimal Latin hypercube designs based on the MaxPro criterion. Parameters in the algorithm may need to be properly tuned to achieve global convergence. Please refer to Joseph, Gul and Ba (2015) for details.
}

\value{
The value returned from the function is a list containing the following components:
\item{Design}{Design matrix}
\item{temp0}{Initial temperature}
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

obj<-MaxProLHD(n = 10, p = 4) 
obj$Design


}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
