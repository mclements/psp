\name{rPR}
\alias{rPR}
\alias{rTPR}
\alias{rFPR}
\title{
  Regression modelling of the relative true positive fraction and
  relative false positive fraction from a paired screen-positive design.
}
\description{
  Regression modelling of the relative true positive fraction and
  relative false positive fraction from a paired screen-positive design.
}
\usage{
rPR(formula, data, test1, test0, TP)
rTPR(formula, data, test1, test0)
rFPR(formula, data, test1, test0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    LHS is the outcome, RHS is the regression equation in terms of the
    covariates in \code{data} and a derived variable \code{test} from
    the two tests.
}
  \item{data}{
    data-frame of the input data
}
  \item{test1}{
    expression for the new test
}
  \item{test0}{
    expression for the reference test
}
  \item{TP}{
    logical to show whether use rTPF (TP=TRUE) or rFPF (TP=FALSE)
}
}
\details{
  Represents the reference test as \code{test} in the regression formula.
}
\value{
  returns a geepack::geeglm object.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
\dontrun{

## Prostate cancer examples from Pepe and Alonzo (2001):

require(foreign)
ex <- read.dta("http://research.fhcrc.org/content/dam/stripe/diagnostic-biomarkers-statistical-center/files/psa_dre_v2.dta")

summary(rTPR(I(d=="yes") ~ test + I(race=="black") + test:I(race=="black"),
             data = ex,
             test1=I(psa=="pos"),
             test0=I(dre=="pos")))
## reduced model
summary(rTPR(I(d=="yes") ~ test + test:I(race=="black"),
             data = ex,
             test1=I(psa=="pos"),
             test0=I(dre=="pos")))

summary(rFPR(I(d=="yes") ~ test + I(race=="black") + test:I(race=="black"),
             data = ex,
             test1=I(psa=="pos"),
             test0=I(dre=="pos")))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
