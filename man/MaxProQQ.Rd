\name{MaxProQQ}
\alias{MaxProQQ}

\title{
Maximum Projection Designs with Quantitative and Qualitative Factors
}
\description{
Generate the maximum projection (MaxPro) design with quantitative and qualitative (QQ) factors. The quantitative factors can be continuous factors and/or discrete numeric factors, and the qualitative factors can be nominal factors and/or ordinal factors. If all factors are continuous, please use \code{\link{MaxProLHD}} directly. 

To use this function, first convert all the ordinal factors into discrete numeric factors through the scoring method (see, e.g., Wu and Hamada 2009, Section 14.10). That is, if the ordinal factor has levels, "poor", "fair" and "good", then depending on the nature of the classification, the experimenter can choose some discrete numeric levels such as (1,2,3) or (1,4,5) to represent the three ordinal levels. Also, all columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1]

The users need to provide their desirable sub-design matrix for the nominal factors as part of the initial design input to this function. This function only optimizes the sub-design matrix of the continuous and discrete numeric factors as well as their row orders in joining with the pre-determined sub-design matrix of the nominal factors. 
}
\usage{
MaxProQQ(InitialDesign, p_nom=0, temp0=0, nstarts=1, itermax=400, total_iter=1e+06)
}

\arguments{
  \item{InitialDesign}{
The initial design matrix where each row is an experimental run and each column is a factor. The rightmost p_nom columns correspond to the p_nom nominal factors, and the columns on the left are for continuous factors and discrete numeric factors (including ordinal factors). It is recommended to use a random Latin hypercube design as the initial design for continuous factors, and all columns of the continuous and discrete numeric factors should be standardized into the unit range of [0,1]. This function only optimizes the element permutations within the columns for continuous and discrete numeric factors. Elements in the last p_nom columns are pre-determined by the user and fixed in this function.  
}
  \item{p_nom}{
Optional, default is ``0''. The number of nominal factors
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
This function utilizes a version of the simulated annealing algorithm to efficiently generate the optimal design for continuous, nominal, discrete numeric and ordinal types of factors, where the design columns for nominal factors are pre-determined by the users and all the other columns are optimized based on the MaxProQQ criterion for the combined design. Parameters in the algorithm may need to be properly tuned to achieve global convergence. Please refer to Joseph, Gul and Ba (2018) for details.
}

\value{
The value returned from the function is a list containing the following components:
\item{Design}{Design matrix}
\item{temp0}{Initial temperature}
\item{measure}{The MaxProQQ criterion measure}
\item{time_rec}{Time to complete the search}
\item{ntotal}{The total number of iterations}
}
\references{
Joseph, V. R., Gul, E., and Ba, S. (2015) "Maximum Projection Designs for Computer Experiments," \emph{Biometrika}, 102, 371-380.

Joseph, V. R., Gul, E., and Ba, S. (2018) "Designing Computer Experiments with Multiple Types of Factors: The MaxPro Approach," \emph{Journal of Quality Technology}.

Wu, C. F. J., and Hamada, M. (2009), \emph{Experiments: Planning, Analysis, and Parameter Design Optimization, 2nd Edition}, New York: Wiley.
}
\author{
Shan Ba <shanbatr@gmail.com> and V. Roshan Joseph <roshan@isye.gatech.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MaxProRunOrder}}, \code{\link{MaxProAugment}}, \code{\link{MaxProLHD}}
}

\examples{
#18-run design
#3 continuous factors, 1 discrete numeric factor (6 levels), 2 nominal factors (3 levels)

#Generate a random Latin hypercube design (each factor standardized into [0,1]) 
#as the initial design for continuous factors
rand_design_part1=apply(matrix(rep(seq(from=0,to=1,length=18),3),ncol=3),2,sample)

#Generate a random initial design (standardized into [0,1]) for discrete numeric factors
rand_design_part2=sample(rep(seq(from=0,to=1,length=6),each=3))

#Construct an optimal design for the two nominal factors
OA_matrix=cbind(rep(1:3,each=6),rep(1:3,6))

#Initial design matrix
InitialDesign=cbind(rand_design_part1,rand_design_part2,OA_matrix)

#Optimize the design based on MaxProQQ criterion
obj=MaxProQQ(InitialDesign, p_nom=2) 
obj$Design


}

\keyword{ Design of Experiments }
\keyword{ Computer Experiments }
