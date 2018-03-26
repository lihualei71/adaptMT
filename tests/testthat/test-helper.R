## context("Formula completion")

## test_that("complete_formula completes the package name (splines) for simple formulas", {
##     tempfun <- function(formula){
##         formula <- complete_pkg(formula)
##         as.character(as.formula(formula))
##     }
##     res <- c("~", "splines::ns(x, df = 8)")
##     formula1 <- ~ ns(x, df = 8)
##     formula2 <- ~ns(x, df = 8)
##     formula3 <- "~ ns(x, df = 8)"
##     formula4 <- "~ns(x, df = 8)"
##     formula5 <- "ns(x, df = 8)"
##     formula6 <- ~ splines::ns(x, df = 8)
##     formula7 <- ~splines::ns(x, df = 8)
##     formula8 <- "~ splines::ns(x, df = 8)"
##     formula9 <- "~splines::ns(x, df = 8)"
##     formula10 <- "ns(x, df = 8)"

##     expect_equal(tempfun(formula1), res)
##     expect_equal(tempfun(formula2), res)
##     expect_equal(tempfun(formula3), res)
##     expect_equal(tempfun(formula4), res)
##     expect_equal(tempfun(formula5), res)
##     expect_equal(tempfun(formula6), res)
##     expect_equal(tempfun(formula7), res)
##     expect_equal(tempfun(formula8), res)
##     expect_equal(tempfun(formula9), res)
##     expect_equal(tempfun(formula10), res)
## })

## test_that("complete_formula completes the package name (splines) for complex formulas", {
##     res1 <- " ~ splines::ns(x, df = 8)"
##     res2 <- " ~  splines::ns(x, df = 8)"    
##     formula1 <- ~ ns(x, df = 8) + splines::ns(y, df = 8)

##     expect_equal(complete_formula(formula1, ""), res1)
##     expect_equal(complete_formula(formula2, ""), res1)
##     expect_equal(complete_formula(formula3, ""), res2)
##     expect_equal(complete_formula(formula4, ""), res1)
##     expect_equal(complete_formula(formula5, ""), res1)
##     expect_equal(complete_formula(formula6, ""), res1)
##     expect_equal(complete_formula(formula7, ""), res1)
##     expect_equal(complete_formula(formula8, ""), res2)
##     expect_equal(complete_formula(formula9, ""), res1)
##     expect_equal(complete_formula(formula10, ""), res1)
## })

## test_that("complete_formula completes the package name mgcv", {
##     res1 <- " ~ splines::ns(x, df = 8)"
##     res2 <- " ~  splines::ns(x, df = 8)"    
##     formula1 <- ~ ns(x, df = 8)
##     formula2 <- ~ns(x, df = 8)
##     formula3 <- "~ ns(x, df = 8)"
##     formula4 <- "~ns(x, df = 8)"
##     formula5 <- "ns(x, df = 8)"
##     formula6 <- ~ splines::ns(x, df = 8)
##     formula7 <- ~splines::ns(x, df = 8)
##     formula8 <- "~ splines::ns(x, df = 8)"
##     formula9 <- "~splines::ns(x, df = 8)"
##     formula10 <- "ns(x, df = 8)"    

##     expect_equal(complete_formula(formula1, ""), res1)
##     expect_equal(complete_formula(formula2, ""), res1)
##     expect_equal(complete_formula(formula3, ""), res2)
##     expect_equal(complete_formula(formula4, ""), res1)
##     expect_equal(complete_formula(formula5, ""), res1)
##     expect_equal(complete_formula(formula6, ""), res1)
##     expect_equal(complete_formula(formula7, ""), res1)
##     expect_equal(complete_formula(formula8, ""), res2)
##     expect_equal(complete_formula(formula9, ""), res1)
##     expect_equal(complete_formula(formula10, ""), res1)
## })
