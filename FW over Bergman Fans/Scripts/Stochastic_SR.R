#### Load and Functions ####

library(ggpubr)
library(latex2exp)
source("Scripts/Functions.R")
load("Data/Matroid_Info.rda")
boundary_cases <- list()

#### Testing Hausdorff Distances ####

load("Data/RandomK_1000_25.rda")
trials <- length(S_x)

trials_per_samp <- 10
data2 <- matrix(0, trials*trials_per_samp, ncol=4)
r <- 0
for (i in 1:trials){
  k0 <- K_x[i] # Which cone (index) is our center point from?
  preImage_Set <- H_NS[which(H_K == k0)] # What is the pre-image set for this cone?
  K <- K_NS[[k0]] # The actual cone based on the index.
  SU <- S_x[[i]] # Sample from B(M).
  FWU_inK <- fwu_inK[[i]] # The point in FW(S) and K furthest from K's boundary.
  tw_min_i <- tw_min[i] # The distance from the boundary of K. 
  for (j in 1:trials_per_samp){
    r <- r + 1
    # Generate uniform noise, and then scale so that it's ~99% of the threshold.
    Z <- matrix(runif(length(SU), min=-1, max=1), nrow=nrow(SU))
    Z <- Z / (sum(diff(Rfast::rowMinsMaxs(Z))) / (1/2*tw_min_i)) *0.99
    if (sum(diff(Rfast::rowMinsMaxs(Z))) %>=% (1/2*tw_min_i)) warning("Too much noise.")
    if ((sum(diff(Rfast::rowMinsMaxs(Z))) / (1/2*tw_min_i)) %<<% 0.99) warning("Not enough noise.")
    SZ <- SU + Z # Add noise
    FW_SZ <- tr.fw.mcf(SZ) # Compute FW set of noisy data.
    FW_SZ_V <- -t(FW_SZ$K_star) # K_star given as min-plus tropical vertices.
    # Project the original optimal FW vertex onto FW(S+Z).
    FWU_inK_proj <- proj.tconv(FWU_inK, FW_SZ_V)
    # Check to ensure the projection actually maps to FW(S+Z).
    if (any(tr.fw.grad(SZ, FWU_inK_proj) %!=% 0)) stop("Projection not in FW(S+Z)")
    # How far did the projection move (in tropical distance).
    dist_moved <- diff(min_max(FWU_inK-FWU_inK_proj))
    if (dist_moved %>>% tw_min_i) stop("Bad Hausdorff") # Checking Hausdorff
    data2[r,] <- c(i, j, tw_min_i, dist_moved)
    # Projection onto FW(S+Z) not in B(M), so project (again) now onto B(M).
    FWU_inK_projU <- proj.tconv(FWU_inK_proj, V_matrix)
    # Test cone membership. Did we stay in K?
    if (which_cone(FWU_inK_projU, K_NS) != k0) stop("Safety Radius Violated")
    print(paste0(i,"/", length(S_x), " : ", j,"/",trials_per_samp))
  }
}

HausDT <- data.table(data2)
colnames(HausDT) <- c("Sample","Trial","Radius","Distance")
save(HausDT, file = "Results/RandomK_1000_25_Haus.rda")

load("Results/RandomK_1000_25_Haus.rda")
HausPlot <- ggplot(HausDT[Distance %>>% 0], aes(x=Radius, y=Distance))+
  geom_point(pch=21, alpha=0.5) + scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color='red', lty=2, lwd=1) +
  ylab(TeX(r"($d_{tr}(\tilde{w}, FW(S+\epsilon))$)")) + 
  xlab(TeX(r"($\tilde{w}_{min}$)")) +
  theme_minimal()

png("Figs/Haus_Plot.png", width = 300, height = 200); HausPlot; dev.off()

#### Stochastic Safety Radius Empirical Computations ####

load("Data/RandomK_1000_25.rda")
trials <- length(S_x)
trials_per_samp <- 1
sig_facs <- c(1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 3, 4, 5)

data <- matrix(0, nrow=length(sig_facs)*trials*trials_per_samp, ncol=5)
r <- 0
for (sig in sig_facs){
  for (i in 1:trials){
    ki <- K_x[i] # Cone index for center of sample.
    K_code <- cone_codes[ki] # Unique code for K in Nested Set Fan.
    Ki <- K_NS[[ki]] # Actual cone.
    preImage_Set <- H_NS[which(H_K == ki)] # Pre-image set for K.
    SU <- S_x[[i]] # The sample itself.
    FWU_inK <- fwu_inK[[i]] # The point in FW(S) and K furthest from K's boundary.
    tw_min_i <- tw_min[i] # Distance to boundary of K.
    for (j in 1:trials_per_samp){
      r <- r + 1
      # Add Gaussian noise, scaled by the sigma-factor and distance to K's boundary.
      Z <- matrix(rnorm(length(SU), mean=0, sd=sig*tw_min_i), nrow=nrow(SU))
      SZ <- SU + Z # Add noise
      # Can't use the new method (LP) because Theorem 13 doesn't apply.
      # Use method outlined in Algorithm 2.
      in_K <- does_fw_int_K(SZ, V_matrix, C_matrix, ki, Ki, K_code, FWU_inK, preImage_Set)
      # What is the ratio of noise to threshold value?
      Z_w_ratio <- sum(diff(Rfast::rowMinsMaxs(Z))) / (1/2*tw_min_i)
      data[r,] <- c(sig, i, j, Z_w_ratio, in_K)
      print(paste0("Sig: ", sig, " : ", i,"/", length(S_x), " : ", j,"/",trials_per_samp))
    }
  }
}

StochDT <- data.table(data)
colnames(StochDT) <- c("Sigma","Sample","Trial","NoiseRatio","InCone")
#save(StochDT, file = "Results/RandomK_1000_25_10_Stoch.rda")

load("Results/RandomK_1000_25_10_Stoch.rda")
StochDT[, .(avg = mean(InCone)), by = list(Sigma)] # Mean by sigma factor.

p1 <- ggplot(StochDT, aes(x=NoiseRatio, y=InCone)) +
  geom_jitter(height=0.025, width=0, alpha=0.5, pch=21) +
  stat_smooth(method = "glm", method.args = list(family = binomial), se=TRUE, color="red")+
  geom_vline(xintercept = 1, lty=2, alpha=0.5) +
  xlab(TeX(r"($\sum_{j=1}^n \|\epsilon_j \|_{\tr} / (\tilde{w}_{min}/2)$)")) +
  ylab(TeX(r"($\exists x \in FW(S+\epsilon) \,|\, \pi_M(x)\in relint(K)$)")) +
  theme_minimal()
  
p2 <- ggplot(StochDT, aes(x=NoiseRatio, y=InCone)) +
  geom_jitter(height=0.025, width=0, alpha=0.5, pch=21) +
  stat_smooth(method = "glm", method.args = list(family = binomial), se=TRUE, color="red") +
  scale_x_log10() +
  geom_vline(xintercept = 1, lty=2, alpha=0.5) +
  xlab(TeX(r"($\sum_{j=1}^n \|\epsilon_j \|_{\tr} / (\tilde{w}_{min}/2)$)")) +
  ylab(TeX(r"($\exists x \in FW(S+\epsilon) \,|\, \pi_M(x)\in relint(K)$)")) +
  theme_minimal()

Safety_Radius_Plot <- ggarrange(p1, p2, nrow=1)
png("Figs/Safety_Radius_Plot.png", width = 700, height = 300); Safety_Radius_Plot; dev.off()

#### Stochastic SR Lower Bound ####

q <- 8
Z <- matrix(rnorm(1000000*q, mean=0, sd=1), ncol=q)
Z_tr <- as.numeric(diff(rowMinsMaxs(Z)))
Vq <- var(Z_tr)
Eq <- mean(Z_tr)

w_min <- 1
n <- 25
sigma <- seq(from=0, to=w_min/(2*n*Eq), length.out=10000)
p_lb <- 1-(Vq/n)/(w_min/(2*n*sigma)-Eq)^2
kdata <- data.table(Sigma=sigma, Pr_LB = p_lb)
ggplot(kdata, aes(x=Sigma, y=Pr_LB))+
  geom_line() + geom_vline(xintercept=w_min/(2*n*Eq), lty=2, lwd=1, col='red') +
  ylim(0,1)
