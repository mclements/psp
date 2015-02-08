psp.univariate <- function(disease, data, test1, test0) {
    disease <- eval(substitute(disease), data)
    test1 <- eval(substitute(test1), data)
    test0 <- eval(substitute(test0), data)
    a <- sum( disease &  test1 &  test0)
    b <- sum( disease &  test1 & !test0)
    c <- sum( disease & !test1 &  test0)
    e <- sum(!disease &  test1 &  test0)
    f <- sum(!disease &  test1 & !test0)
    g <- sum(!disease & !test1 &  test0)
    rTPR <- (a+b)/(a+c)
    varlogrTPR <- (b+c)/(a+b)/(a+c)
    rFPR <- (e+f)/(e+g)
    varlogrFPR <- (f+g)/(e+f)/(e+g)
    selogrTPR <- sqrt(varlogrTPR)
    selogrFPR <- sqrt(varlogrFPR)
    cirTPR <- rTPR*exp(c(qnorm(0.025)*selogrTPR,
                         qnorm(0.975)*selogrTPR))
    cirFPR <- rFPR*exp(c(qnorm(0.025)*selogrFPR,
                         qnorm(0.975)*selogrFPR))
    test.rTPR <- stats:::mcnemar.test(matrix(c(0,b,c,0),2))
    test.rFPR <- stats:::mcnemar.test(matrix(c(0,f,g,0),2))
    structure(list(a=a,b=b,c=c,e=e,f=f,g=g,
                   rTPR=rTPR,
                   rFPR=rFPR,
                   varlogrTPR=varlogrTPR,
                   varlogrFPR=varlogrFPR,
                   cirTPR=cirTPR,
                   cirFPR=cirFPR,
                   p.rTPR=test.rTPR$p.value,
                   p.rFPR=test.rFPR$p.value),
              class="psp.univ")
}
print.psp.univ <- function(obj, ...) {
    cat(sprintf("\nEstimated rTPR=%f (95%% CI: %f, %f; p=%f)\n",
                obj$rTPR,
                obj$cirTPR[1],
                obj$cirTPR[2],
                obj$p.rTPR))
    cat(sprintf("\nEstimated rFPR=%f (95%% CI: %f, %f; p=%f)\n",
                obj$rFPR,
                obj$cirFPR[1],
                obj$cirFPR[2],
                obj$p.rFPR))
    invisible(obj)
}
