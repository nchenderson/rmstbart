\name{rmstbart}
\alias{rmstbart}
\title{
RMST BART for survival analysis
}
\description{
RMST-BART performs inference on the restricted mean survival time (RMST)
function.For a survival time \eqn{T}, covariate vector \eqn{x}, and positive
number \eqn{\tau}, the RMST function is defined as
\eqn{E[ \min(T, \tau)| x] = \mu(x)}

RMST-BART models the unknown RMST
function \eqn{\mu(x)} using a Bayesian \dQuote{sum-of-trees} model.\cr

\eqn{\mu(x)} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{\mu(x)}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a small amount to the overall fit.
}
\usage{
rmstbart(times, delta, x.train, Gweights, x.test=matrix(0.0,0,0),
    tau=NULL, sgrid=NULL, sigma.mu=NULL,
    transformation="identity",
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=TRUE,
    cont=FALSE, rm.const=TRUE,
    kappa=2.0, power=2.0, base=.95, sigmaf=NA,
    ntree=200L, numcut=100L,
    ndpost=1000L, nskip=100L, keepevery=1L,
    nkeeptrain=ndpost, nkeeptest=ndpost,
    nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
    printevery=100L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{times}{
     The time of event or time of right-censoring.
   }
   \item{delta}{
     The event indicator: 1 is an event while 0 is censored.
   }
   \item{x.train}{ Explanatory variables for training (in sample)
    data.\cr May be a matrix or a data frame, with (as usual) rows
    corresponding to observations and columns to variables.\cr If a
    variable is a factor in a data frame, it is replaced with dummies.
    Note that \eqn{q} dummies are created if \eqn{q>2} and one dummy
    created if \eqn{q=2} where \eqn{q} is the number of levels of the
    factor.  \code{rmstbart} will generate draws of \eqn{f(x)} for each
    \eqn{x} which is a row of \code{x.train}.  }

   \item{x.test}{ Explanatory variables for test (out of sample)
   data. Should have same structure as \code{x.train}.
   \code{rmstbart} will generate draws of \eqn{f(x)} for each \eqn{x} which
   is a row of \code{x.test}.  }

   \item{transformation}{
      Function to transform survival times. Can be equal to \code{"identity"}
      or \code{"log"}
   }


   \item{sparse}{Whether to perform variable selection based on a
     sparse Dirichlet prior rather than simply uniform; see Linero 2016.}
   \item{theta}{Set \eqn{theta} parameter; zero means random.}
   \item{omega}{Set \eqn{omega} parameter; zero means random.}
   \item{a}{Sparse parameter for \eqn{Beta(a, b)} prior:
     \eqn{0.5<=a<=1} where lower values inducing more sparsity.}
   \item{b}{Sparse parameter for \eqn{Beta(a, b)} prior; typically,
     \eqn{b=1}.}
   \item{rho}{Sparse parameter: typically \eqn{rho=p} where \eqn{p} is the
     number of covariates under consideration.}
   \item{augment}{Whether data augmentation is to be performed in sparse
     variable selection.}

   \item{xinfo}{ You can provide the cutpoints to BART or let BART
     choose them for you.  To provide them, use the \code{xinfo}
     argument to specify a list (matrix) where the items (rows) are the
     covariates and the contents of the items (columns) are the
     cutpoints.  }

   \item{usequants}{ If \code{usequants=FALSE}, then the
    cutpoints in \code{xinfo} are generated uniformly; otherwise,
    if \code{TRUE}, uniform quantiles are used for the cutpoints. }

   \item{rm.const}{ Whether or not to remove constant variables.}

   \item{sigest}{Preliminary estimate of residual}

   \item{kappa}{ For numeric \eqn{y}, \code{k} is the number of prior
   standard deviations \eqn{E(Y|x) = f(x)} is away from +/-0.5.  The
   response, code{y.train}, is internally scaled to range from -0.5 to
   0.5.  For binary \eqn{y}, \code{kappa} is the number of prior standard
   deviations \eqn{f(x)} is away from +/-3.  The bigger \code{kappa} is, the more
   conservative the fitting will be.  }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{sigmaf}{
    The SD of \eqn{f}.  Not used if \eqn{y} is binary.
   }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{numcut}{ The number of possible values of \eqn{c} (see
   \code{usequants}).  If a single number if given, this is used for all
   variables.  Otherwise a vector with length equal to
   \code{ncol(x.train)} is required, where the \eqn{i^{th}}{i^th}
   element gives the number of \eqn{c} used for the \eqn{i^{th}}{i^th}
   variable in \code{x.train}.  If usequants is false, numcut equally
   spaced cutoffs are used covering the range of values in the
   corresponding column of \code{x.train}.  If \code{usequants} is true, then
   \eqn{min(numcut, the number of unique values in the corresponding
   columns of x.train - 1)} values are used.  }

   \item{ndpost}{
   The number of posterior draws returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.
   }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \code{rmstbart} returns an object of class \code{rmstbart}.
  This is a list with the following named components

   \item{yhat.train}{
   A matrix with ndpost rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with ndpost rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
\keyword{survival}
\keyword{tree}
\keyword{nonparametric}
\keyword{regression}
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
