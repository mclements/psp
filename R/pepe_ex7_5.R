require(geepack)
lhs <- function (formula) 
    if (length(formula) == 3) formula[[2]] else NULL
rhs <- function (formula) 
    if (length(formula) == 3) formula[[3]] else formula[[2]]
`lhs<-` <- function (formula, value) 
{
    if (length(formula) == 2) 
        as.formula(as.call(c(formula[[1]], value, formula[[2]])))
    else {
        newformula <- formula
        newformula[[2]] <- value
        newformula
    }
}
`rhs<-` <- function (formula, value) 
{
    newformula <- formula
    newformula[[length(formula)]] <- value
    newformula
}
rPR <- function(formula, data, test1, test0, TP) {
    disease <- eval(lhs(formula), data)
    test1 <- eval(substitute(test1), data)
    test0 <- eval(substitute(test0), data)
    data <- transform(data, .id = 1:nrow(data))
    newdata <- rbind(transform(data,
                               test=1,
                               tp=ifelse(disease==0 | is.na(disease),0,1)*test1, 
                               fp=ifelse(disease==1 | is.na(disease),0,1)*test1),
                     transform(data,
                               test=0,
                               tp=ifelse(disease==0 | is.na(disease),0,1)*test0,
                               fp=ifelse(disease==1 | is.na(disease),0,1)*test0))
    newdata <- newdata[order(newdata$.id), ]
    lhs(formula) <- if (TP) quote(tp) else quote(fp)
    geeglm(formula,
           data = newdata,
           family = binomial(link="log"),
           id = .id)
}
rTPR <- function(...) rPR(..., TP = TRUE)
rFPR <- function(...) rPR(..., TP = FALSE)
    
require(foreign)
ex <- read.dta("~/work/psa_dre_v2.dta")

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

