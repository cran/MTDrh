\name{MTDrh}
\alias{MTDrh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Mass Transportation Distance Rank Histogram
}
\description{
Constructs a mass transportation distance rank histogram to assess the reliability of probabilistic scenarios using observations for a set of instances [1].
}
\usage{
MTDrh(scenarios, observation, prob = NULL, debias = FALSE, transformation = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scenarios}{ A dataset that contains scenarios. It should be a 3 dimensional array:
%%     ~~Describe \code{scenarios} here~~
(dimension of each scenario)x(number of scenarios per instance)x(number of instances)
}
  \item{observation}{A dataset that contains observations. The dimension of each observation and the number of instances should match the dimension and number of instances of the scenarios. It should be a matrix:

%%     ~~Describe \code{observation} here~~
(dimension of each observation)x(number of instances)
}
  \item{prob}{ A dataset that contains the probability of each scenario for each instance.
  If prob is not given, the default that the scenarios have equal probabilities.
%%     ~~Describe \code{prob} here~~
It should be a matrix:
  (number of scenarios)x(number of instances)
}
  \item{debias}{ If debias=TRUE, the data are debiased for each instance [1].
%%     ~~Describe \code{debias} here~~
}
  \item{transformation}{If transformation=TRUE, the data are transformed with Mahalanobis transformation for each instance [1].
%%     ~~Describe \code{transformation} here~~
}
}

\value{
%%  ~Describe the value returned
Returns an array of mass transportation ranks and a histogram plot.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
[1] D. Sari, Y. Lee, S. Ryan, D. Woodruff. Statistical metrics for assessing the quality of
wind power scenarios for stochastic unit commitment. Wind Energy 
19, 873-893 (2016)
doi:10.1002/we.1872
}
\author{ Didem Sari, Sarah M. Ryan
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
#Generate 1000 instances of 10 scenarios and observation with dimension 8 
#from the same normal distribution.

scen <- array(rnorm(8*10*1000,0,1),dim=c(8,10,1000)) 
obs <- array(rnorm(8*1000,0,1),dim=c(8,1000))
ranks <- MTDrh(scen,obs,prob=NULL,debias=FALSE,transformation=FALSE)

#Generate 1000 instances of 27 scenarios and observation with dimension 8 
#from AR(1) processes. The marginal distributions of the scenarios and observation
#are the same but the autocorrelation levels are different. The Mahalanobis
#transformation is applied. See Figure 8 [1].

scen <- array(arima.sim(list(order=c(1,0,0),ar=0.10),n=8*27*1000,sd=1),dim=c(8,27,1000))
obs <- array(arima.sim(list(order=c(1,0,0),ar=0.90),n=8*1000,sd=0.45),dim=c(8,1000))
ranks<-MTDrh(scen,obs,prob=NULL,debias=FALSE,transformation=TRUE)
hist(ranks, breaks=c(0:28),xlab="bin",ylab="frequency",col="gray",main="MTD rh")

#Generate 1000 instances of 27 scenarios that have heterogeneous autocorrelation
#levels and corresponding observations with autocorrelation different
#from the scenarios.
#The marginal standard deviations of scenarios and observation match.  See Figure 9 [1]

scen1 <- array(arima.sim(list(order=c(1,0,0),ar=0.10),n=8*10*1000,sd=1),dim=c(8,10,1000))
scen2 <- array(arima.sim(list(order=c(1,0,0),ar=0.80),n=8*17*1000,sd=0.64),dim=c(8,17,1000))
scen <- array(NA,dim=c(8,27,1000))
scen[,1:10,]<-scen1
scen[,11:27,]<-scen2
obs <- array(arima.sim(list(order=c(1,0,0),ar=0.50),n=8*1000,sd=0.86),dim=c(8,1000))
ranks<-MTDrh(scen,obs,prob=NULL,debias=FALSE,transformation=TRUE)
hist(ranks, breaks=c(0:28),xlab="bin",ylab="frequency",col="gray",main="MTD rh")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mass Transportation Distance rank histogram }% use one of  RShowDoc("KEYWORDS")
