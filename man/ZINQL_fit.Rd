\name{ZINQL_fit}
\alias{ZINQL_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ZINQL_fit
}
\description{
A Zero-Inflated Quantile Approach for Differential Abundance Analysis of Longitudinal Microbiome Data. The test combines both logistic mixed-effects model and longitudinal quantile regression model. The test is applied to individual taxon.
}
\usage{
ZINQL_fit(y, formula, formula.logistics = NA,
  meta, C, taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
  method = "both", n.positive.cut = 5,
  seed = 2024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The taxa read count. It should be numeric.}
  \item{formula}{The full model for ZINQ-L. By default, this is the formula for both logistics mixed-effects regression and longitudinal quantile regression. For example, y ~ X + Y + (1|ID), where y is zero-inflated.}
  \item{formula.logistics}{The model of logistics mixed-effects regression. By default the user does not need to input, and this will be the same as \code{formula}.}
  \item{meta}{A data.frame including all covariates of interest. All of them should be numeric. The user should transfer factor variables to binary or numeric values.}
  \item{C}{The name of clinical variable of interest. The final P-value will tell whether the abundance is associated with this variable.}
  \item{taus}{A grid of quantile levels, e.g., 0.5 for the median, 0.75 for the 3rd quartile; default is c(0.1, 0.25, 0.5, 0.75, 0.9).}
  \item{method}{Combination method, "MinP" for MinP combination, "Cauchy" for Truncated Cauchy combination, or 'Both' which will return results for both methods. The default is "Both".}
  \item{n.positive.cut}{Minimum number of positive values for y. If the condition is not satisfied, the test will not be conducted and NA P-values will be returned.}
  \item{seed}{A seed for any steps with randomness. The default is 2024.}
}
\details{
\itemize{
\item ZINQ-L requires at least \code{n.positive.cut} of positive \code{y}. If this condition is not satisfied, the test will not be conducted and NA P-values will be returned.
\item By default, the user does not need to input \code{formula.logistics}. This means that the formula of logistics mixed-effects regression is consistent with the longitudinal quantile regression. If the user wants the formula of them to be different, the user can input a different formula for logistics mixed-effects regression, and the \code{formula} will be the model for longitudinal quantile regression.
\item ZINQ-L will automatically drop covariates with unique values. ZINQ-L will also drop variables that have unique values on positive part of \code{y} for longitudinal quantile regression.
}
}
\value{
The function returns a list:
\itemize{
  \item Final_P_value - This is the Final P-value of ZINQ-L.
  \item model - This is a String variable suggesting the status of the test. If the test is successful conducted, this will be either 'ZINQL_MinP', 'ZINQL_Cauchy' or 'Both'. If the test fails, this will be String 'None'.
  \item Intermediate_P_value - These are intermediate P-values that are combined to form the final P-value by either Truncated Cauchy or MinP approach.
}
}
\references{
references
}
\author{
The package is developed by Shuai Li, Runzhe Li, Wodan Ling and Ni Zhao.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#generate example data
y = c(rep(0,5),1,2,5,6,7)
meta = data.frame(
  ID=c('1','1','1','2','2','2','3','4','3','4'),
  gender=c(0,0,0,1,1,1,0,1,0,1),
  x=c(0,0,0,1,1,1,1,0,1,0),
  age=c(1,2,3,2,3,4,1,2,2,3)
)

#run the function
result = ZINQL_fit(y=y, formula=y~x+age+gender+(1|ID),
  C='x', taus=c(0.1, 0.25, 0.5, 0.75, 0.9), seed=2024,
  meta=meta, method='Both')
print( result )
}
