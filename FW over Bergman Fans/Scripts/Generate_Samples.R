#### Loading and Functions ####

source("Scripts/Functions.R")
load("Data/Matroid_Info.rda")
boundary_cases <- list()

# Function for Algorithm 1.
S_int_K <- function(x0, k0, n, dispersion=3){
  K <- K_NS[[k0]]
  X0 <- Rfast::rep_row(x0, n)
  w_min_x0 <- w_min_C(x0, C_matrix)
  # Generate and apply noise.
  Z <- matrix(rnorm(length(X0), mean=0, sd=dispersion*w_min_x0), nrow=n)
  S <- round(X0 + Z, 6) # Rounding helps avoid issues with FW computations.
  # Project the points back into U.
  SU <- t(apply(S, 1, function(x) proj.tconv(x, V_matrix)))
  if (any(!(apply(SU, 1, is_M_ultrametric, C_matrix=C_matrix)))) stop("BAD PROJECTION")
  # Compute maximal w_min for any x in the intersection of FW(SU) with K.
  x_opt <- max_w_min(K, SU, C_matrix)
  while (is.null(x_opt)) {
    Z <- matrix(rnorm(length(X0), mean=0, sd=dispersion*w_min_x0), nrow=n)
    S <- round(X0 + Z, 6) # Rounding helps avoid issues with FW computations.
    SU <- t(apply(S, 1, function(x) proj.tconv(x, V_matrix)))
    if (any(!(apply(SU, 1, is_M_ultrametric, C_matrix=C_matrix)))) stop("BAD PROJECTION")
    x_opt <- max_w_min(K, SU, C_matrix)
  }
  # Some potentially useful stats.
  SU_K <- t(apply(SU, 1, function(x) recover_K(x)$Nested)) # Distribution of cones for sample.
  SU_inK <- length(which(SU_K[,k0] %==% 1)) # How many were in K?
  SU_meanK <- mean(SU_K[,k0])
  tw_min <- w_min_C(x_opt, C_matrix)
  FWU_inK_maxWmin <- x_opt
  return(list(S=SU, S_inK=SU_inK, S_meanK=SU_meanK, 
    FWU_inK=FWU_inK_maxWmin, tw_min=tw_min))
}

#### Example 1 ####

nrow(B_matrix) # Number of bases.
nrow(C_matrix) # Number of circuits.
nrow(F_matrix) # Number of (proper) flats.
length(K_BF) # Number of maximal cones in Bergman Fan.
length(K_NS) # Number of maximal cones in Nested Set Fan.
length(K_CF) # Number of maximal cones in Cyclic Bergman Fan.

w <- c(5,5,1,2,2,0,3,3)
K_BF[which_cone(w, K_BF)] # Which cone of Bergman Fan?
K_S1 <- K_NS[[which_cone(w, K_NS)]] # Which cone of Nested Set Fan?
lambda <- c(2,2,2,3)
mu <- 5 # Show construction of x using nested structure.
mu + Rfast::colsums(t(Rfast::eachrow(t(K_S1), lambda))) 

x <- c(5,5,1,2,2,0,6,6)
which_cone(x, K_BF) == which_cone(w, K_BF)
which_cone(x, K_NS) == which_cone(w, K_NS)
which_cone(x, K_CF) == which_cone(w, K_CF)
# Which rays are shared or distinguished by the cyclic Bergman cones?
fintersect(data.table(K_CF[[which_cone(w,K_CF)]]), data.table(K_CF[[which_cone(x,K_CF)]]))
fsetdiff(data.table(K_CF[[which_cone(w,K_CF)]]), data.table(K_CF[[which_cone(x,K_CF)]]))
F_df[is_cyclic & !connected]

#### Example 2 ####

eps <- 3/4*c(1,-1,-1,1,1,-1,-1,-1)
w_prime <- w + eps
w_prime_proj <- proj.tconv(w_prime, V_matrix)
which_cone(w_prime_proj, K_NS) == which_cone(w, K_NS)
which_cone(w_prime_proj, K_CF) == which_cone(w, K_CF)

#### Example 3 ####

S <- rbind(v1=w, v2=c(5,5,1,2,2,0,4,4), v3=x)
eps <- matrix(c(5,5.01,5.01,5,4.99,5,1.01,1,0.99,2,1.99,2.01,2.01,2.01,2,-0.01,0,0,2.99,4,5,2.99,4,4.99), 3, 8)
Se <- S+eps
FW_Se_LP <- tr.fw.mcf(Se)
FW_Se <- norm.i(FW_Se_LP$K_star) # Normalize each entry to canonical coordinates.
unique(FW_Se) == Se[2,]

is_M_ultrametric(Se[2,], C_matrix) # Is it in the Bergman fan?
v2_prime_proj <- proj.tconv(Se[2,], V_matrix) # So we need to project it.
v2_prime_proj == S[2,] # Is it equal?

#### Sample Center Points around Fan ####

num_samples <- 1000
sample_size <- 25

K_x <- sample(seq_len(length(K_NS)), num_samples, replace = TRUE)
X <- t(sapply(K_x, function(k) sample_cone(K_NS[[k]])))

S_x <- list()
S_n_inK <- c()
S_avg_inK <- c()
fwu_inK <- list()
tw_min <- c()
for (i in 1:num_samples){
  xi <- as.numeric(X[i,])
  ki <- K_x[[i]]
  find_sample <- S_int_K(x0=xi, k0=ki, n=sample_size, dispersion=3)
  S_x[[length(S_x)+1]] <- find_sample$S
  S_n_inK <- c(S_n_inK, find_sample$S_inK)
  S_avg_inK <- c(S_avg_inK, find_sample$S_meanK)
  fwu_inK[[length(fwu_inK)+1]] <- find_sample$FWU_inK
  tw_min <- c(tw_min, find_sample$tw_min)
  print(paste0(i,"/", num_samples, " Samples Generated"))
}
print(paste0("Avg S in K: ", round(mean(S_avg_inK),3), "; Range of tw_min: ", paste0(round(min_max(tw_min),4), collapse=" - "),3))

save(K_x, X=X, S_x, S_n_inK, S_avg_inK, fwu_inK, tw_min, file = "Data/RandomK_1000_25.rda")
