#### Load and Functions ####

library(ggpubr)
library(latex2exp)
source("Scripts/Functions.R")
load("Data/M(Running)_Info.rda")

#### Testing Hausdorff Distance for Uniform Samples (Appendix) ####

trials <- 1000
n_vec <- c(2,4,8,16)
q_vec <- c(2,4,6,8)
HD_mat <- matrix(0, trials*length(n_vec)*length(q_vec), 4)
r <- 0
for (q in q_vec){
  for (n in n_vec){
    t <- 0
    while (t < trials){
      S <- matrix(runif(n*q), n, q) # Generate initial sample S.
      Z <- matrix(rnorm(n*q), n, q) # Generate noise Z.
      FW_S <- tr.fw.lp(S, full_set = TRUE)$FW.set # FW(S).
      FW_Z <- tr.fw.lp(S+Z, full_set = TRUE)$FW.set # FW(S+Z)
      if (!is.null(FW_S)){
        if (any(!apply(FW_S, 1, function(x) all(tr.fw.grad(S, x) %==% 0)))) FW_S <- NULL
      } else next
      if (!is.null(FW_Z)) {
        if (any(!apply(FW_Z, 1, function(x) all(tr.fw.grad(S+Z, x) %==% 0)))) FW_Z <- NULL
      } else next
      FW_Str <- -t(tr.fw.mcf(S)$K_star) # Max-plus tropical vertices of FW(S)
      FW_Ztr <- -t(tr.fw.mcf(S+Z)$K_star) # Max-plus tropical vertices of FW(S+Z)
      if (!is.null(FW_Str)){
        if (any(!apply(FW_Str, 1, function(x) all(tr.fw.grad(S, x) %==% 0)))) FW_Str <- NULL
      } else next
      if (!is.null(FW_Ztr)) {
        if (any(!apply(FW_Ztr, 1, function(x) all(tr.fw.grad(S+Z, x) %==% 0)))) FW_Ztr <- NULL
      } else next
      # If there were no errors in computations.
      if (!is.null(FW_S) & !is.null(FW_Z) & !is.null(FW_Str) & ! is.null(FW_Ztr)){
        t <- t + 1
        r <- r + 1
        FW_S_proj <- t(apply(FW_S, 1, function(x) proj.tconv(x, FW_Ztr))) # Project onto FW(S+Z). This is a min tropical dist mapping.
        FW_Z_proj <- t(apply(FW_Z, 1, function(x) proj.tconv(x, FW_Str))) # Project onto FW(S)
        d_SZ <- as.numeric(diff(Rfast::rowMinsMaxs(FW_S_proj - FW_S))) # Compute tropical distance to projection.
        d_ZS <- as.numeric(diff(Rfast::rowMinsMaxs(FW_Z_proj - FW_Z))) # Ditto.
        hd_dist <- max(c(d_SZ, d_ZS)) # Take the maximum over all distances.
        z_tr <- sum(diff(Rfast::rowMinsMaxs(Z))) # Compute the sum of tropical norms of noise.
        HD_mat[r,] <- c(q, n, hd_dist, z_tr)
        print(paste0("q = ", q, " : n = ", n, " : ", t, " / ", trials))
      }
    }
  }
}
HD_dt <- data.table(HD_mat)
colnames(HD_dt) <- c("q","n","d","z")
HD_dt[, scaled_shift := d/z]
#save(HD_dt, file = "Results/U(0,1)_1000_Hausdorff_App.rda")

load("Results/U(0,1)_1000_Hausdorff_App.rda")
min_max(HD_dt$scaled_shift)
HD_dt[which.max(HD_dt$scaled_shift)]

(HD_demo <- ggplot(HD_dt, aes(x = scaled_shift, y = after_stat(density))) +
  geom_histogram(binwidth = 0.05) + 
  facet_grid(n~q, labeller = label_both) +
  xlab("Scaled Shift") +
  ylab("Density") +
  theme_bw()) +
  theme(axis.text = element_text(size=10),
    axis.title = element_text(size=14),
    strip.text = element_text(size=12))

png("Figs/Haus_Demo.png", width = 600, height = 400); HD_demo; dev.off()

#### Testing Hausdorff Distances for M-Ultra Samples of M ####

load("Data/M(Running)_S1000_N25.rda")
num_samples <- length(S_x)

trials_per_samp <- 10
SU_HD_mat <- matrix(0, num_samples*trials_per_samp, ncol=5)
r <- 0
for (i in 1:num_samples){
  SU_i <- S_x[[i]] # Sample from B(M).
  w_i <- w_mat[i,] # The point in FW(S) and K furthest from K's boundary.
  w_min_i <- w_min_vec[i] # The distance from the boundary of K. 
  j <- 0
  while (j < trials_per_samp){
    # Generate uniform noise, and then scale so that it's at the threshold.
    Z <- matrix(runif(length(SU_i), min=-1, max=1), nrow=nrow(SU_i))
    Z <- Z / (sum(diff(Rfast::rowMinsMaxs(Z))) / (w_min_i/2)) # Scale the noise.
    SZ <- SU_i + Z # Add noise
    FW_SZ <- tr.fw.mcf(SZ) # Compute FW set of noisy data.
    if (is.null(FW_SZ)) next # In case dual LP is degenerate.
    r <- r + 1
    j <- j + 1
    P <- -t(FW_SZ$K_star) # K_star given as min-plus tropical vertices.
    # Project the original optimal FW vertex onto FW(S+Z).
    piP_w <- proj.tconv(w_i, P)
    # Check to ensure the projection actually maps to FW(S+Z).
    if (any(tr.fw.grad(SZ, piP_w) %!=% 0)) stop("Projection not in FW(S+Z)")
    # How far did the projection move (in tropical distance).
    dtr.w.piP_w <- diff(min_max(w_i - piP_w))
    if (dtr.w.piP_w %>>% w_min_i) stop("Bad Hausdorff") # Checking Hausdorff
    # Projection onto FW(S+Z) not in B(M), so project (again) now onto B(M).
    piU_piP_w <- proj.tconv(piP_w, V_matrix)
    # Test cone membership. Did we stay in K?
    dtr.w.piU_piP_w <- diff(min_max(w_i - piU_piP_w))
    if (dtr.w.piU_piP_w %>=% w_min_i) stop("Safety Radius Violated")
    SU_HD_mat[r,] <- c(i, j, w_min_i, dtr.w.piP_w, dtr.w.piU_piP_w)
    print(paste0(i,"/", length(S_x), " : ", j,"/",trials_per_samp))
  }
}

SU_HD_DT <- data.table(SU_HD_mat)
colnames(SU_HD_DT) <- c("Sample","Trial","w_min","d_H","Shift")
#save(SU_HD_DT, file = "Results/M(Running)_S1000_N25_R10_HD.rda")

load("Results/M(Running)_S1000_N25_R10_HD.rda")

SU_HD_DT[, `:=` (wH_Ratio = d_H/w_min, wS_Ratio = Shift/w_min)]
min_max(SU_HD_DT$wH_Ratio)
min_max(SU_HD_DT$wS_Ratio)
#load("Results/RandomK_1000_25_Haus.rda")
(SU_HausPlot <- ggplot(SU_HD_DT[d_H %>>% 0], aes(x=w_min, y=d_H))+
  geom_point(pch=21, alpha=0.5) + 
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color='red', lty=2, lwd=1) +
  ylab(TeX(r"($d_{tr}(\tilde{w}, FW(S+\epsilon))$)")) + 
  xlab(TeX(r"($\tilde{w}_{min}$)")) +
  theme_bw() +
  theme(axis.text = element_text(size=10),
    axis.title = element_text(size=16),
    strip.text = element_text(size=12)))

png("Figs/Haus_Plot_K5.png", width = 500, height = 350); SU_HausPlot; dev.off()



#### Stochastic Safety Radius Empirical Computations ####

load("Data/M(Running)_S1000_N25.rda")
num_samples <- length(S_x)
trials_per_samp <- 1
sig_facs <- c(1/32, 1/16, 1/8, 1/4, 1/2, 1, 2, 3, 4, 5)

SSR_mat <- matrix(0, nrow=length(sig_facs)*num_samples*trials_per_samp, ncol=5)
r <- 0
for (sig in sig_facs){
  for (i in 1:num_samples){
    ki <- K_x[i] # Cone index for center of sample.
    K_code <- cone_codes[ki] # Unique code for K in Nested Set Fan.
    Ki <- K_NS[[ki]] # Actual cone.
    preImage_Set <- H_NS[which(H_K == ki)] # Pre-image set for K.
    SU <- S_x[[i]] # The sample itself.
    w_i <- w_mat[i,] # The point in FW(S) and K furthest from K's boundary.
    w_min_i <- w_min_vec[i] # Distance to boundary of K.
    for (j in 1:trials_per_samp){
      r <- r + 1
      # Add Gaussian noise, scaled by the sigma-factor and distance to K's boundary.
      Z <- matrix(rnorm(length(SU), mean=0, sd=sig*w_min_i), nrow=nrow(SU))
      SZ <- SU + Z # Add noise
      # Use method outlined in Algorithm 1.
      in_K <- does_fw_int_K(SZ, V_matrix, C_matrix, ki, Ki, K_code, w_i, preImage_Set)
      # What is the ratio of noise to threshold value?
      Z_w_ratio <- sum(diff(Rfast::rowMinsMaxs(Z))) / (w_min_i/2)
      SSR_mat[r,] <- c(sig, i, j, Z_w_ratio, in_K)
      print(paste0("Sig: ", sig, " : ", i,"/", length(S_x), " : ", j,"/",trials_per_samp))
    }
  }
}

StochDT <- data.table(SSR_mat)
colnames(StochDT) <- c("Sigma","Sample","Trial","NoiseRatio","InCone")
#save(StochDT, file = "Results/M(Running)_S1000_N25_V10_R1_SSR.rda")

load("Results/M(Running)_S1000_N25_V10_R1_SSR.rda")
StochDT[, .(avg = mean(InCone)), by = list(Sigma)] # Mean by sigma factor.

(p1 <- ggplot(StochDT, aes(x=NoiseRatio, y=InCone)) +
  geom_jitter(height=0.025, width=0, alpha=0.5, pch=21) +
  stat_smooth(method = "glm", method.args = list(family = binomial), se=TRUE, color="red")+
  geom_vline(xintercept = 1, lty=2, alpha=0.5) +
  xlab("Noise Ratio") +
  ylab("Trial Outcome") + 
  #xlab(TeX(r"($\sum_{j=1}^n \|\epsilon_j \|_{\tr} / (\tilde{w}_{min}/2)$)")) +
  #ylab(TeX(r"($\exists x \in FW(S+\epsilon) \,|\, \pi_M(x)\in relint(K)$)")) +
  theme_bw() +
  theme(axis.text = element_text(size=10),
    axis.title = element_text(size=14)))
  
(p2 <- ggplot(StochDT, aes(x=NoiseRatio, y=InCone)) +
  geom_jitter(height=0.025, width=0, alpha=0.5, pch=21) +
  stat_smooth(method = "glm", method.args = list(family = binomial), se=TRUE, color="red") +
  scale_x_log10() +
  geom_vline(xintercept = 1, lty=2, alpha=0.5) +
  xlab("Noise Ratio") +
  ylab("Trial Outcome") + 
  #xlab(TeX(r"($\sum_{j=1}^n \|\epsilon_j \|_{\tr} / (\tilde{w}_{min}/2)$)")) +
  #ylab(TeX(r"($\exists x \in FW(S+\epsilon) \,|\, \pi_M(x)\in relint(K)$)")) +
  theme_bw() +
  theme(axis.text = element_text(size=10),
    axis.title = element_text(size=14)))

Safety_Radius_Plot <- ggarrange(p1, p2, nrow=1)
Safety_Radius_Plot <- ggarrange(SU_HausPlot, p1, nrow=1)
png("Figs/Safety_Radius_Plot.png", width = 800, height = 300); Safety_Radius_Plot; dev.off()

#### All Matroids on [q=8] ####

# Read in the matroid bases from the output of SageMath.
matroids.txt <- readLines("Data/all_8_element_matroids.txt")
all_matroids <- list()
bases <- list()
for (line in matroids.txt){
  if (substr(line, 1, 7) == "Matroid") {
    if (length(bases) > 0) all_matroids[[length(all_matroids)+1]] <- bases
    bases <- list()
  }
  if (substr(line, 1, 9) == "frozenset"){
    b <- regmatches(line, regexpr("\\{[^}]+\\}", line))
    if (length(b) > 0){
      b <- gsub("[\\{\\}]", "", b)
      b <- as.numeric(unlist(strsplit(b, ",")))
      # Don't include any of rank 0, 1, or 8.
      if (length(b) > 1 & length(b) < 8) bases[[length(bases)+1]] <- b+1
    }
  }
}
loopless <- sapply(all_matroids, function(x){
  x <- sort(unique(unlist(x)))
  if (length(x) < 8) return(FALSE)
  return(all(x == seq(8)))})
coloops <- sapply(all_matroids, function(x){
  y <- lapply(x, char_vect, E=seq(8))
  z <- do.call(rbind, y)
  any(colsums(z) == length(y))})

loopless_matroids <- all_matroids[which(loopless & !coloops)]

num_samples <- 10
sample_size <- 25
trials_per_samp <- 1
sig_facs <- c(1/32, 1/8, 1, 3, 5)

SSR_mat <- matrix(0, length(all_matroids)*num_samples*length(sig_facs)*trials_per_samp, 6)
m <- 0
r <- 0
for (matroid in loopless_matroids[213:length(loopless_matroids)]) {
  m <- m+1
  # Matroid Info
  M <- matroid_info(matroid)
  E_set <- M$E_set
  B_matrix <- M$B_matrix
  C_matrix <- M$C_matrix
  if (!M_connected(C_matrix)) next
  F_matrix <- M$F_matrix
  V_matrix <- M$V_matrix
  F_df <- M$F_df
  K_BF <- M$K_BF
  K_NS <- M$K_N
  if (any(sapply(K_NS, function(x) is.null(nrow(x))))) next
  K_CF <- M$K_CF
  cone_codes <- M$cone_codes
  H_K <- M$H_K
  H_NS <- M$H_NS
  # Generate Samples
  K_x <- sample(seq_len(length(K_NS)), num_samples, replace = TRUE)
  X <- t(sapply(K_x, function(k) sample_cone(K_NS[[k]])))
  S_x <- list()
  w_mat <- matrix(0, nrow=num_samples, ncol=length(E_set))
  w_min_vec <- c()
  for (i in 1:num_samples){
    v0 <- as.numeric(X[i,])
    K <- K_NS[[K_x[i]]]
    Sample.from.U_M <- Generate.Sample.from.U_M(v0, K, sample_size, sigma_s=3, in_U=TRUE)
    S_x[[length(S_x)+1]] <- Sample.from.U_M$S
    w_mat[i,] <- Sample.from.U_M$w
    w_min_vec <- c(w_min_vec, Sample.from.U_M$w_min)
    print(paste0(m,"/",length(loopless_matroids)," : ", i,"/", num_samples, " Samples Generated"))
  }
  # Stochastic Safety Radius
  for (sig in sig_facs){
    for (i in 1:num_samples){
      ki <- K_x[i] # Cone index for center of sample.
      K_code <- cone_codes[ki] # Unique code for K in Nested Set Fan.
      Ki <- K_NS[[ki]] # Actual cone.
      preImage_Set <- H_NS[which(H_K == ki)] # Pre-image set for K.
      SU <- S_x[[i]] # The sample itself.
      w_i <- w_mat[i,] # The point in FW(S) and K furthest from K's boundary.
      w_min_i <- w_min_vec[i] # Distance to boundary of K.
      for (j in 1:trials_per_samp){
        r <- r + 1
        # Add Gaussian noise, scaled by the sigma-factor and distance to K's boundary.
        Z <- matrix(rnorm(length(SU), mean=0, sd=sig*w_min_i), nrow=nrow(SU))
        SZ <- SU + Z # Add noise
        # Use method outlined in Algorithm 1.
        in_K <- does_fw_int_K(SZ, V_matrix, C_matrix, ki, Ki, K_code, w_i, preImage_Set)
        # What is the ratio of noise to threshold value?
        Z_w_ratio <- sum(diff(Rfast::rowMinsMaxs(Z))) / (w_min_i/2)
        SSR_mat[r,] <- c(m, sig, i, j, Z_w_ratio, in_K)
        print(paste0("Sig: ", sig, " : ", i,"/", length(S_x), " : ", j,"/",trials_per_samp))
      }
    }
  }
}

SSR_final <- SSR_mat[SSR_mat[,1] > 0,]
length(unique(SSR_final[,1]))

allMatroid_DT <- data.table(SSR_final)
colnames(allMatroid_DT) <- c("Matroid","Sigma","Sample","Trial","NoiseRatio","InCone")
#save(allMatroid_DT, file = "Results/M(All_8)_S10_N25_V5_R1.rda")

allMatroid_DT[, .(avg = mean(InCone)), by = list(Sigma)] # Mean by sigma factor.

(allMat_plot <- ggplot(allMatroid_DT, aes(x=NoiseRatio, y=InCone, group=Matroid)) +
    geom_jitter(height=0.025, width=0, alpha=0.5, pch=21) +
    #stat_smooth(method = "glm", method.args = list(family = binomial), 
    #  se=FALSE, color="red", lwd=0.5, alpha=0.01) +
    geom_line(stat="smooth", method = "glm", method.args = list(family=binomial),
      se=FALSE, color="red", lwd=0.5, alpha=0.1) +
    geom_line(data=allMatroid_DT[Matroid==1134], mapping=aes(x=NoiseRatio, y=InCone, group=Matroid), 
      stat="smooth", method = "glm", method.args = list(family = binomial), 
        se=FALSE, color="blue", lwd=1, alpha=1) +
    scale_x_log10() +
    geom_vline(xintercept = 1, lty=2, alpha=0.5) +
    xlab("Noise Ratio") +
    ylab("Trial Outcome") + 
    #xlab(TeX(r"($\sum_{j=1}^n \|\epsilon_j \|_{\tr} / (\tilde{w}_{min}/2)$)")) +
    #ylab(TeX(r"($\exists x \in FW(S+\epsilon) \,|\, \pi_M(x)\in relint(K)$)")) +
    theme_bw() +
    theme(axis.text = element_text(size=10),
      axis.title = element_text(size=14)))

png("Figs/all_matroids.png", width = 600, height = 350); allMat_plot; dev.off()
