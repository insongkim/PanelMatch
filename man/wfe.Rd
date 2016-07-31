\name{wfe}

\alias{wfe}
\title{Fitting the Weighted Fixed Effects Model for Causal Inference}

\description{
  \code{wfe} is used to fit weighted fixed effects model for causal
  inference. \code{wfe} also derives the regression weights for
  different causal quantity of interest.
}

\usage{
wfe(formula, data, treat = "treat.name",
    unit.index, time.index = NULL, method = "unit",
    qoi = c("ate", "att") , estimator = NULL, C.it = NULL,
    hetero.se = TRUE, auto.se = TRUE,
    White = TRUE, White.alpha = 0.05,
    verbose = TRUE, unbiased.se = FALSE, unweighted = FALSE,
    store.wdm = FALSE,
    tol = sqrt(.Machine$double.eps))
}



\arguments{

  \item{formula}{a symbolic description of the model to be fitted. The
    formula should not include dummmies for fixed effects. The details
    of model specifications are given under `Details'.
  }
  \item{data}{data frame containing the variables in the model.
  }
  \item{treat}{a character string indicating the name of treatment
    variable used in the models. The treatment should be binary
    indicator (integer with 0 for the control group and 1 for the
    treatment group).
  }
  \item{unit.index}{a character string indicating the name of unit
    variable used in the models. The index of unit should be factor.
  }
  \item{time.index}{a character string indicating the name of time
    variable used in the models. The index of time should be factor.
  }
  \item{method}{method for weighted fixed effects regression, either
    \code{unit} for unit fixed effects; \code{time} for time fixed
    effects. The default is \code{unit}. 
  }

  \item{qoi}{one of \code{"ate"} or \code{"att"}. The default is
    \code{"ate"}.
  }

  \item{estimator}{an optional character string indicating the
  estimating method. One of \code{"fd"} or \code{"did"}. The default is
  \code{NULL}.
  }
  
  \item{C.it}{an optional non-negative numeric vector specifying
    relative weights for each unit of analysis.
  }
  
  \item{hetero.se}{a logical value indicating whether heteroskedasticity
    across units is allowed in calculating standard errors. The default
    is \code{TRUE}.
  }

  \item{auto.se}{a logical value indicating whether arbitrary
    autocorrelation is allowed in calculating standard errors. The
    default is \code{TRUE}.
  }
 
  \item{White}{a logical value indicating whether White misspecification
    statistics should be calculated. The default is \code{TRUE}.
  }
  
  \item{White.alpha}{level of functional specification test. See White
    (1980) and Imai and Kim (2012). The default is \code{0.05}.
  }

  \item{verbose}{logical. If \code{TRUE}, helpful messages along with
    a progress report of the weight calculation are printed
    on the screen. The default is \code{TRUE}.
  }

  \item{unbiased.se}{logical. If \code{TRUE}, bias-asjusted
    heteroskedasticity-robust standard errors are used. See Stock and
    Watson (2008). Should be used only for balanced panel. The default
    is \code{FALSE}.
  }
  
  \item{unweighted}{logical. If \code{TRUE}, standard unweighted fixed
    effects model is estimated. The default is \code{FALSE}. Note: users
    do not need to specify \code{qoi} when unweighted=\code{TRUE}.}

  \item{store.wdm}{logical. If \code{TRUE}, weighted demeaned
    dataframe will be stored. The default is \code{FALSE}.}

  \item{tol}{a relative tolerance to detect zero singular values for
    generalized inverse. The default is sqrt(.Machine$double.eps)
  }
  

  
}


\value{

 \code{wfe} returns an object of class "wfe", a list that contains the
 components listed below.


 The function \code{summary} (i.e., \code{summary.wfe}) can be used to
 obtain a table of the results.

 \item{coefficients}{a named vector of coefficients}
 \item{residuals}{the residuals, that is respons minus fitted values}
 \item{df}{the degree of freedom}
 \item{W}{weight matrix calculated from the model. Row and column
 indices can be found from \code{units}, \code{times}.}
 \item{call}{the matched call}
 \item{causal}{causal quantity of interest}
 \item{estimator}{the estimating method}
 \item{units}{a dataframe containing unit names used for \code{W}}
 \item{times}{a dataframe containing time names used for \code{W}}
 \item{method}{call of the method used}
 \item{vcov}{the variance covariance matrix}
 \item{White.alpha}{the alpha level for White specification test}
 \item{White.pvalue}{the p-value for White specification test}
 \item{White.stat}{the White statistics}
 \item{X}{the design matrix}
 \item{Y}{the response vector}
 \item{X.wdm}{the demeaned design matrix}
 \item{Y.wdm}{the demeaned response vector}
 \item{mf}{the model frame}
}

\details{

  To fit the weighted unit (time) fixed effects model, use the syntax
  for the formula, \code{y ~ x1 + x2}, where \code{y} is a dependent
  variable and \code{x1} and \code{x2} are unit (time) varying
  covariates.

  \code{wfe} calculates weights based on different underlying causal
  quantity of interest: Average Treatment Effect (\code{qoi = "ate"}) or
  Average Treatment Effect for the Treated (\code{qoi = "att"}).

  One can further set estimating methods: First-Difference
  (\code{estimator ="fd"}) or Difference-in-differences (\code{estimator
    = "did"}). For the two-way fixed effects model, set \code{estimator
    = "did"}
  
  To specify different ex-ante weights for each unit of analysis, use
  non-negative weights \code{C.it}. For instance, using the survey
  weights for \code{C.it} enables the estimation fo the average
  treatement effect for the target population.

  An object of class "wfe" contains vectors of unique unit(time) names and
  unique unit(time) indices. 

}

\references{
  Imai, Kosuke and In Song Kim. (2012) ``On the Use of Linear Fixed
  Effects Regression Models for Causal Inference.'' Technical Report,
  Department of Politics, Princeton University.  available at
  \url{http://imai.princeton.edu/research/FEmatch.html}

  Stock, James and Mark Watson. (2008) ``Heteroskedasticity-Robust
  Standard Errors for Fixed Effect Panel Data Regression'' Econometrica,
  76, 1.
  
  White, Halbert. (1980) `Using Least Squares to Approximate Unknown
  Regression Functions.''  International Economic Review, 21, 1,
  149--170.

}
\seealso{
  \code{pwfe} for fitting weighted fixed effects models with propensity
  score weighting
}
\examples{
### NOTE: this example illustrates the use of wfe function with randomly
### generated panel data with arbitrary number of units and time.

## generate panel data with number of units = N, number of time = Time
N <- 10 # number of distinct units
Time <- 15 # number of distinct time

## treatment effect
beta <- 1

## generate treatment variable
treat <- matrix(rbinom(N*Time, size = 1, 0.25), ncol = N)
## make sure at least one observation is treated for each unit
while ((sum(apply(treat, 2, mean) == 0) > 0) | (sum(apply(treat, 2, mean) == 1) > 0) |
       (sum(apply(treat, 1, mean) == 0) > 0) | (sum(apply(treat, 1, mean) == 1) > 0)) {
  treat <- matrix(rbinom(N*Time, size = 1, 0.25), ncol = N)
}
treat.vec <- c(treat)

## unit fixed effects
alphai <- rnorm(N, mean = apply(treat, 2, mean))

## geneate two random covariates
x1 <- matrix(rnorm(N*Time, 0.5,1), ncol=N)
x2 <- matrix(rbeta(N*Time, 5,1), ncol=N)
x1.vec <- c(x1)
x2.vec <- c(x2)
## generate outcome variable
y <- matrix(NA, ncol = N, nrow = Time)
for (i in 1:N) {
    y[, i] <- alphai[i] + treat[, i] + x1[,i] + x2[,i] + rnorm(Time)
}
y.vec <- c(y)

## generate unit and time index
unit.index <- rep(1:N, each = Time)
time.index <- rep(1:Time, N)

Data.str <- as.data.frame(cbind(y.vec, treat.vec, unit.index, x1.vec, x2.vec))
colnames(Data.str) <- c("y", "tr", "strata.id", "x1", "x2")

Data.obs <- as.data.frame(cbind(y.vec, treat.vec, unit.index, time.index, x1.vec, x2.vec))
colnames(Data.obs) <- c("y", "tr", "unit", "time", "x1", "x2")

############################################################
# Example 1: Stratified Randomized Experiments
############################################################

## run the weighted fixed effect regression with strata fixed effect.
## Note: the quantity of interest is Average Treatment Effect ("ate")
## and the standard errors allow heteroskedasticity and arbitrary
## autocorrelation.


### Average Treatment Effect
mod.ate <- wfe(y~ tr+x1+x2, data = Data.str, treat = "tr",
               unit.index = "strata.id", method = "unit",
               qoi = "ate", hetero.se=TRUE, auto.se=TRUE)
## summarize the results
summary(mod.ate)

### Average Treatment Effect for the Treated
mod.att <- wfe(y~ tr+x1+x2, data = Data.str, treat = "tr",
               unit.index = "strata.id", method = "unit",
               qoi = "att", hetero.se=TRUE, auto.se=TRUE)
## summarize the results
summary(mod.att)


############################################################
# Example 2: Observational Studies with Unit Fixed-effects
############################################################

## run the weighted fixed effect regression with unit fixed effect.
## Note: the quantity of interest is Average Treatment Effect ("ate")
## and the standard errors allow heteroskedasticity and arbitrary
## autocorrelation.

mod.obs <- wfe(y~ tr+x1+x2, data = Data.obs, treat = "tr",
               unit.index = "unit", time.index = "time", method = "unit",
               qoi = "ate", hetero.se=TRUE, auto.se=TRUE,
               White = TRUE, White.alpha = 0.05)

## summarize the results
summary(mod.obs)

## extracting weigths
summary(mod.obs)$W



###################################################################
# Example 3: Observational Studies with differences-in-differences
###################################################################

## run difference-in-differences estimator.
## Note: the quantity of interest is Average Treatment Effect ("ate")
## and the standard errors allow heteroskedasticity and arbitrary
## autocorrelation.

mod.did <- wfe(y~ tr+x1+x2, data = Data.obs, treat = "tr",
               unit.index = "unit", time.index = "time", method = "unit",
               qoi = "ate", estimator ="did", hetero.se=TRUE, auto.se=TRUE,
               White = TRUE, White.alpha = 0.05, verbose = TRUE)

## summarize the results
summary(mod.did)

## extracting weigths
summary(mod.did)$W

}


\author{Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
  and In Song Kim, Princeton University, \email{insong@princeton.edu}}

\keyword{regression}
