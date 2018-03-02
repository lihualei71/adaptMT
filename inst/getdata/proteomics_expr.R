source("other_methods.R")
source("AdaPT.R")
library("IHWpaper")
library("IHW")

proteomics_file <- system.file(
    "extdata/real_data",
    "science_signaling.csv",
    package = "IHWpaper"
    )

proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)

proteomics_df$pvalue <- rank(
    proteomics_df$p1,
    ties.method="first"
    ) * proteomics_df$p1 / nrow(proteomics_df) 

pvals <- proteomics_df$pvalue
x <- log(proteomics_df$X..peptides)

reorder <- rev(order(x))
x <- x[reorder]
pvals <- pvals[reorder]
x <- data.frame(x = x)

alphalist <- seq(0.01, 0.3, by=0.01) # target FDR level
tau <- 0.5; eps <- 0.1 # parameters for SABHA
thr <- 0.5 # parameter for Storey-BH
thr1 <- 0.1; thr2 <- 0.5 # parameters for adaptive SeqStep
ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
set.seed(1)
n <- length(pvals)

num_alpha <- length(alphalist)
max_alpha <- max(alphalist)
# gather results
NumRej <- matrix(0, nrow = 13, ncol = num_alpha)

# methods: 1 SeqStep, 2 HingeExp, 3 ForwardStop, 4 Adaptive SeqStep, 5 BH, 6 Storey-BH, 7 Barber-Candes, 8 SABHA (step), 9 SABHA (ordered), 10 IHW, 11 IHW (oracle), 12 IF (oracle), 13 AdaPT

qhat_step <- Solve_q_step(pvals, tau, eps)
qhat_ordered <- Solve_q_ordered_simple(pvals, tau, eps, ADMM_params)
for(i in 1:num_alpha){
    NumRej[1,i] <- SeqStep(pvals, alpha = alphalist[i], C = 2)
    NumRej[2,i] <- HingeExp(pvals, alpha = alphalist[i])
    NumRej[3,i] <- ForwardStop(pvals, alpha = alphalist[i])
    NumRej[4,i] <- length(Adaptive_SeqStep_method(pvals, alphalist[i], alphalist[i], tau))
    NumRej[5,i] <- length(BH_method(pvals, alphalist[i]))
    NumRej[6,i] <- length(Storey_method(pvals, alphalist[i], thr))
    NumRej[7,i] <- length(BC_method(pvals, alphalist[i]))
    NumRej[8,i] <- length(SABHA_method(pvals, qhat_step, alphalist[i], tau))
    NumRej[9,i] <- length(SABHA_method(pvals, qhat_ordered, alphalist[i], tau))
    NumRej[10,i] <- rejections(ihw(pvals, x[, 1], alphalist[i]))
    NumRej[11,i] <- ihw.oracle(pvals, 1:n, alphalist[i])
    NumRej[12,i] <- IF.oracle(pvals, x[, 1], alphalist[i])
}

print("AdaPT starts!")

pi.formulas <- paste0("ns(x, df = ", 6:10, ")")
mu.formulas <- paste0("ns(x, df = ", 6:10, ")")
formulas <- expand.grid(pi.formulas, mu.formulas)

res.AdaPT <-
    AdaPT.glm(x, pvals,
              dist = beta.family(),
              pi.formulas = formulas[, 1],
              mu.formulas = formulas[, 2])

NumRej[13,] <- res.AdaPT$num_rej[1:num_alpha]

result <- list(NumRej = NumRej,
               res.AdaPT = res.AdaPT)

save(file = "data/proteomics_res.RData", result)
