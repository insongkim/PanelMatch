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
    dyad1.index = NULL, dyad2.index = NULL,
    qoi = "ate", estimator = NULL, C.it = NULL,
    hetero.se = TRUE, auto.se = TRUE, df.adjustment = TRUE,  
    dyad.se = FALSE,
    White = TRUE, White.alpha = 0.05,
    verbose = TRUE, unbiased.se = FALSE, unweighted = FALSE,
    store.wdm = FALSE, maxdev.did = NULL,
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
    effects. The default is \code{unit}. For two-way weighted fixed
    effects regression models, set method to the default value
    \code{unit}.}

  \item{dyad1.index}{a character string indicating the variable name of first unit
    of a given dyad. The default is \code{NULL}.
    This is required to calculate robust standard errors with dyadic data.
  }

  \item{dyad2.index}{a character string indicating the variable name of second unit
    of a given dyad. The default is \code{NULL}.
    This is required to calculate robust standard errors with dyadic data.
  }

  \item{qoi}{one of \code{"ate"} or \code{"att"}. The default is
    \code{"ate"}. If set to \code{"att"} in implementing \code{"fd"}
    and \code{"did"} estimators, the comparison of the treated
    observation is restricted to the control observation from the
    previous time period but not with the control observation from the
    next time period.}

  \item{estimator}{an optional character string indicating the
  estimating method. One of \code{"fd"}, \code{"did"}, or
  \code{"Mdid"}. \code{"fd"} is for First-Difference
  Design. \code{"did"} is for multi-period Difference-in-Differences
  design. The default is \code{NULL}. Setting estimator to be
  \code{"Mdid"} implements the Difference-in-Differences design with
  Matching on the pretreatment outcome variables.}
  
  \item{C.it}{an optional non-negative numeric vector specifying
    relative weights for each unit of analysis. If not specified, the
    weights will be calculated based on the estimator and quantity of
    interest. }
  
  \item{hetero.se}{a logical value indicating whether heteroskedasticity
    across units is allowed in calculating standard errors. The default
    is \code{TRUE}.
  }

  \item{auto.se}{a logical value indicating whether arbitrary
    autocorrelation is allowed in calculating standard errors. The
    default is \code{TRUE}.
  }

  \item{dyad.se}{a logical value indicating whether correlations across dyads exist. The
    default is \code{FALSE}. 
  }

  \item{df.adjustment}{a logical value indicating whether a degrees of freedom adjustment 
      for standard errors calculation is performed. The default is \code{TRUE}.
  }

  \item{White}{a logical value indicating whether White misspecification
    statistics should be calculated. The default is \code{TRUE}.
  }
  
  \item{White.alpha}{level of functional specification test. See White
    (1980) and Imai and Kim (2016). The default is \code{0.05}.
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
    effects model is estimated. The default is \code{FALSE}. Note:
    users do not need to specify \code{qoi} when
    unweighted=\code{TRUE}. For standard two-way fixed effects model
    (unit and time), set estimator=\code{"did"} and
    unweighted=\code{"TRUE"}.}

  \item{store.wdm}{logical. If \code{TRUE}, weighted demeaned
    dataframe will be stored. The default is \code{FALSE}.}

  \item{maxdev.did}{an optional positive numeric value specifying the
  maximum deviation in pre-treatment outcome when \code{"Mdid"} is
  implemented. The default is \code{NULL}, which implements
  nearest-neighbor matching.}

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
 \item{W}{a dataframe containing unit and time indices along with the
 weights used for the observation. If method=\code{unit}, integer numbers
 corresponding to the order of input data will be used for generating time index.}
 \item{Num.nonzero}{Number of observations with non-zero weights}
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
 \item{mf}{the model frame where the last column is the weights used for the analysis}
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

\references{Imai, Kosuke and In Song Kim. (2016) ``When Should We Use
  Linear Fixed Effects Regression Models for Causal Inference with
  Longitudinal Data?'' Technical Report, Department of Politics,
  Princeton University.  available at
  \url{http://imai.princeton.edu/research/FEmatch.html}

  Aronow, Peter M., Cyrus Samii, and Valentina A. Assenova (2015) ``Cluster–robust
  Variance Estimation for Dyadic Data." Political Analysis 23, no. 4, 564--577.
	
  Stock, James and Mark Watson. (2008) ``Heteroskedasticity-Robust
  Standard Errors for Fixed Effect Panel Data Regression'' Econometrica,
  76, 1.
  
  White, Halbert. (1980) ``Using Least Squares to Approximate Unknown
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


\dontrun{
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

#########################################################################
# Example 4: DID with Matching on Pre-treatment Outcomes
#########################################################################

## implements matching on pre-treatment outcomes where the maximum
## deviation is specified as 0.5 

mod.Mdid <- wfe(y~ tr+x1+x2, data = Data.obs, treat = "tr",
               unit.index = "unit", time.index = "time", method = "unit",
               qoi = "ate", estimator ="Mdid", hetero.se=TRUE, auto.se=TRUE,
               White = TRUE, White.alpha = 0.05, maxdev.did = 0.5, verbose = TRUE)

## summarize the results
summary(mod.Mdid)

## Note: setting the maximum deviation to infinity (or any value
## bigger than the maximum pair-wise difference in the outcome) will
## return the same result as Example 3.

dev <- 1000+max(Data.obs$y)-min(Data.obs$y)
mod.did2 <- wfe(y~ tr+x1+x2, data = Data.obs, treat = "tr",
               unit.index = "unit", time.index = "time", method = "unit",
               qoi = "ate", estimator ="Mdid", hetero.se=TRUE, auto.se=TRUE,
               White = TRUE, White.alpha = 0.05, maxdev.did = dev, verbose = TRUE)

## summarize the results
summary(mod.did2)
mod.did2$coef[1] ==  mod.did$coef[1]
}

}


\author{In Song Kim, Massachusetts Institute of Technology, \email{insong@mit.edu}		
	and Kosuke Imai, Princeton University, \email{kimai@princeton.edu}	}

\keyword{regression}
