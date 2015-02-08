\name{psp.univariate}
\alias{psp.univariate}
\title{
  Univariate estimation of the relative true positive fraction and
  relative false positive fraction from a paired screen-positive design.
}
\description{
  Univariate estimation of the relative true positive fraction and
  relative false positive fraction from a paired screen-positive design.
  Variance estimation uses 
}
\usage{
psp.univariate(disease, data, test1, test0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disease}{
    expression for the outcome
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
}
\details{
  This uses the variance estimator from Cheng and Macaluso (1997).
}
\value{
  returns a list of class \code{psp.univ} with elements
  \code{rTPR},\code{varlogrTPR},\code{rFPR} and \code{varlogrFPR}. At
  present, the only method available is \code{print.psp.univ}, which
  prints the point estimates and confidence intervals.
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

## Prostate cancer example from Pepe and Alonzo (2001):

require(foreign)
ex <- read.dta("http://research.fhcrc.org/content/dam/stripe/diagnostic-biomarkers-statistical-center/files/psa_dre_v2.dta")

psp.univariate(disease=I(d=="yes"), data=ex,
               test1=I(psa=="pos"),
               test0=I(dre=="pos"))

## We can also do this, albeit less efficiently, using regression models:

summary(fit <- rTPR(I(d=="yes") ~ test,
                    data = ex,
                    test1=I(psa=="pos"),
                    test0=I(dre=="pos")))
summary(fit <- rFPR(I(d=="yes") ~ test,
                    data = ex,
                    test1=I(psa=="pos"),
                    test0=I(dre=="pos")))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
