psp.univariate1 <- function(a,b,c,level=0.95,label="rPR") {
    rPR <- (a+b)/(a+c)
    varlogrPR <- (b+c)/(a+b)/(a+c)
    selogrPR <- sqrt(varlogrPR)
    aa <- (1 - level)/2
    aa <- c(aa, 1 - aa)
    fac <- stats::qnorm(aa)
    ci <- rPR*exp(fac*selogrPR)
    test <- stats:::mcnemar.test(matrix(c(0,b,c,0),2))
    structure(list(a=a,b=b,c=c,
                   statistic=rPR,
                   ci=ci,
                   test=test,
                   level=level,
                   label=label),
              class="psp.univariate1")
}
print.psp.univariate1 <- function(obj, label=NULL, eps=1e-5, ...) {
    if (!is.null(label)) obj$label <- label
    mat <- with(obj, matrix(c(a,b,c,NA),2,byrow=TRUE))
    dimnames(mat) <- list("Test 1"=c(TRUE,FALSE),"Test 2"=c(TRUE,FALSE))
    print.table(mat)
    with(obj, 
         cat(sprintf("\n%s=%f (%i%% CI: %f, %f; p-value: %s)\n",
                     label,
                     statistic,
                     level*100,
                     ci[1],
                     ci[2],
                     format.pval(test$p.value,eps=eps))))
    invisible(obj)
}
psp.univariate <- function(disease, data, test1, test0, level = 0.95) {
    disease <- eval(substitute(disease), data)
    test1 <- eval(substitute(test1), data)
    test0 <- eval(substitute(test0), data)
    a <- sum( disease &  test1 &  test0)
    b <- sum( disease &  test1 & !test0)
    c <- sum( disease & !test1 &  test0)
    e <- sum(!disease &  test1 &  test0)
    f <- sum(!disease &  test1 & !test0)
    g <- sum(!disease & !test1 &  test0)
    rTPR <- psp.univariate1(a,b,c,level,label="rTPR")
    rFPR <- psp.univariate1(e,f,g,level,label="rFPR")
    structure(list(rTPR=rTPR, rFPR=rFPR),
              class="psp.univariate")
}
print.psp.univariate <- function(obj, ...) {
    cat("Disease:\n")
    print(obj$rTPR, ...)
    cat("\nNo disease:\n")
    print(obj$rFPR, ...)
    invisible(obj)
}

if (FALSE) {
    library(psp)
    psp.univariate1(179,264,137,level=0.9)
    psp.univariate1(138,717,978,level=0.9)
    expand.data <- function(data,n)
        data[rep(1:nrow(data), eval(substitute(n), data)), , drop=FALSE]
    d <- expand.data(data.frame(d=c(1,1,1,0,0,0),
                                test1=c(1,1,0,1,1,0),
                                test2=c(1,0,1,1,0,1),
                                n=c(179,264,137,138,717,978)), n)
    print(psp.univariate(d,d,test1,test2,level=0.9), eps=1e-10)
}

