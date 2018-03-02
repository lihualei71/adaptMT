library('dplyr')
source("AdaPT.R")
source("other_methods.R")

#####################################################################################
## run fMRI data example
#####################################################################################

alphalist <- seq(0.01, 0.3, by<-0.01) # target FDR level
tau <- 0.5 # parameters for SABHA
thr <- 0.5 # parameter for Storey-BH
ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr

set.seed(1)
load("fMRI.RData")
labels <- as.factor(output$ROI_labels)
pvals <- output$pvals
n <- length(pvals)

x <- data.frame(output$voxel_coords, output$ROI_labels)
names(x) <- c("x", "y", "z", "labels")
x <- mutate(x, x2 = x^2, y2 = y^2, z2 = z^2, labels = labels)
x <- model.matrix(~ labels + 0 + (x + y + z)^2 + x2 + y2 + z2, data = x)

## x <- data.frame(output$voxel_coords)
## names(x) <- c("x", "y", "z")
## x <- mutate(x, x2 = x^2, y2 = y^2, z2 = z^2)
## x <- model.matrix(~ (x + y + z)^2 + x2 + y2 + z2 + 0, data = x)

num_alpha <- length(alphalist)
max_alpha <- max(alphalist)
# gather results
NumRej <- matrix(0, nrow = 8, ncol = num_alpha)

# methods: 1 BH, 2 Storey-BH, 3 Barber-Candes, 4 SABHA (eps = 0.1), 5 SABHA (eps = 0.2), 6 SABHA (eps = 0.5), 7 IHW, 8 AdaPT
SABHA.factor.root <- 1 / 2 / (1 - tau) * sum(sqrt(summary(labels))) / n

eps1 <- 0.1
qhat_block1 <- Solve_q_block(pvals, tau, eps1, labels, ADMM_params)
SABHA.factor1 <- 1 + SABHA.factor.root / eps1

eps2 <- 0.2
qhat_block2 <- Solve_q_block(pvals, tau, eps2, labels, ADMM_params)
SABHA.factor2 <- 1 + SABHA.factor.root / eps2

eps3 <- 0.5
qhat_block3 <- Solve_q_block(pvals, tau, eps3, labels, ADMM_params)
SABHA.factor3 <- 1 + SABHA.factor.root / eps3

for(i in 1:num_alpha){
    NumRej[1,i] <- length(BH_method(pvals, alphalist[i]))
    NumRej[2,i] <- length(Storey_method(pvals, alphalist[i], thr))
    NumRej[3,i] <- length(BC_method(pvals, alphalist[i]))
    NumRej[4,i] <- length(SABHA_method(pvals, qhat_block1, alphalist[i] / SABHA.factor1, tau))
    NumRej[5,i] <- length(SABHA_method(pvals, qhat_block2, alphalist[i] / SABHA.factor2, tau))
    NumRej[6,i] <- length(SABHA_method(pvals, qhat_block3, alphalist[i] / SABHA.factor3, tau))
    NumRej[7,i] <- rejections(ihw(pvals, labels, alphalist[i]))
}

res.AdaPT <- AdaPT.glmnet(x = x, pvals = pvals,
                          dist = beta.family(),
                          glmnet.args = list(nlambda = 30))
NumRej[8,] <- res.AdaPT$num_rej[1:30]

result <- list(NumRej = NumRej,
               res.AdaPT = res.AdaPT)
save(file = "data/fMRI_res.RData", result)
