\name{CoxSummary}
\alias{CoxSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of Cox proportional hazards regression
}
\description{
Display the outcome of multivariate Cox proportional hazards regression in a dataframe.
}
\usage{
CoxSummary(models, var, summ, styled)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
A list of models. Each model is an output of coxph.
}
\item{var}{
A string with the name of the variable of interest.
}
\item{summ}{
Boolian of whether the output should be summarized.
}
\item{styled}{
Boolian of whether the output should be styled.
}
}

\references{
Therneau T (2023). A Package for Survival Analysis in R. R package version 3.5-5, https://CRAN.R-project.org/package=survival.
}
\author{
Tamas Szili-Torok
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
CoxSummary (models = list(model.1, model.2), var = 'variable', summ = TRUE, styled = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
