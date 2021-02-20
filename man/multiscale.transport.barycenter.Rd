\name{multiscale.transport.barycenter}
\alias{multiscale.transport.barycenter}
\title{Multiscale Transport Barycenter}
\description{
Compute optimal transport barycenter using a multiscale linear programming approach.
}


\usage{

multiscale.transport.barycenter(X.mean, weights.mean = 1, gmra.marginals, weights.marginals , trp.lp, p=2, nType=0, dType=1, scaleMass =TRUE, mean.radius=10^-10, step=0.5, n.iterations=100, ret.sol=TRUE, store.sol = FALSE, store.prefix="./trp-")

}

\arguments{
  \item{X.mean}{Initial mean point set}
  \item{weights.mean}{Vector of weights for each point in X.mean or a single scalar if all points have the same weight}
  \item{gmra.marginals}{ gmra trees of point sets for which to compute the mean }
  \item{weight/marginals}{Vector of weight or isngle scalar for each point set ni gmra.marginals for weighting the individual points. This is not for weighting the marginal point to compute a weighted mean.}
  \item{trp.lp}{Multiacle lp solver to use for solving the optimal tranpsort problems}
  \item{p}{p-th Wasssertstein distance}
  \item{nType}{Neighborhood type, 0=absolut, 1 = relative}
  \item{dType}{Cost function type. 1 = Euclidean metric, 2 = L1 metric, 3 =
    Squared Euclidean Metric, 4 = Gaussian 2-Wasserstein per Node}
  \item{scaleMass}{Scale the masses of each marginal and the mean to sum to one. If false and}
  \item{mean.radius}{ The radius for the mean gmra tree, i.e., the resolution at which the optimal transport to the mean is computed }
  \item{step}{Stepszie of gradient decent at each iteration for refining the mean}
  \item{ret.sol}{Return the solution transport plans form each marginal to the mean}
  \item{store.sol}{Store the solution on harddisk}
  \item{store.prefix}{Store each transport plan in prefix\{number of marginal\}.Rdata file.}
  \item{n.iteratios}{Maximal number of iterations}
}
\value{
A list with named entries trp.plans, gmra.mean, X.mean and costs (total cost of all transports to the mean at each iteration) }

\author{
  Samuel Gerber
}
\seealso{
  multiscale.transport
}
\references{

}
\examples{

#create example data  
phi <- seq(0, 2*pi, length.out=78)
X1<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0)
X1[,1]=X1[,1]*0.5 + 3;
X1[,2]=X1[,2]*2;

#phi <- runif(353)*2*pi
phi <- seq(0, 2*pi, length.out=79)
X2<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0)
X2[,1]=X2[,1]*3;
X2[,2]=X2[,2]*0.5 + 3;

phi <- seq(0, 2*pi, length.out=76)
X3<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0)
X3[,1]=X3[,1]*0.5 - 3
X3[,2]=X3[,2]*2 

phi <- seq(0, 2*pi, length.out=75)
X4<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0)
X4[,1]=X4[,1]*0.5 
X4[,2]=X4[,2]*2 -4

#create multiscale decompsotions
library(gmra)
gmra1 = gmra.create.ipca(X=X1, eps=0, d=2, maxKids=2)
gmra2 = gmra.create.ipca(X=X2, eps=0, d=2, maxKids=2)
gmra3 = gmra.create.ipca(X=X3, eps=0, d=2, maxKids=2)
gmra4 = gmra.create.ipca(X=X4, eps=0, d=2, maxKids=2)

#setup and solve multiscale lp
library(mop)
trp.lp <- multiscale.transport.create.lp(oType=30)
icprop <- multiscale.transport.create.iterated.capacity.propagation.strategy(3, 0)
multiscale.transport.set.propagation.strategy.1(trp.lp, icprop);

X.mean = cbind(runif(423), runif(423) ) 
#X.mean = cbind(runif(1000), runif(1000) )
#X.mean = rbind(X1)
gmra.marginals = list(gmra1, gmra2, gmra3, gmra4)
weight.marginals = list(1,1, 1, 1)

avg = multiscale.transport.barycenter( X.mean, 1, gmra.marginals, weight.marginals, trp.lp, rel.err=0, p=2, dType=1  )

plot( rbind(X1, X2, X3, X4), asp=1)
points(avg$X.mean, col="red")

}
\keyword{multiscale,gmra,neighbor}
