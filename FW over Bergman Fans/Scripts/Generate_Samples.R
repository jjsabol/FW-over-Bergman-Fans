#### Loading and Functions ####

source("Scripts/Functions.R")
load("Data/M(Running)_Info.rda")

# Function to Generate Sample from M-Ultrametric Space
Generate.Sample.from.U_M <- function(v0, K, n, sigma_s=3, in_U=FALSE){
  V0 <- Rfast::rep_row(v0, n)
  w_min_v0 <- w_min_C(v0, C_matrix)
  # Generate and apply noise.
  Z <- matrix(rnorm(length(V0), mean=0, sd=sigma_s*w_min_v0), nrow=n)
  S <- V0 + Z
  # Project the points back into U.
  if (in_U) { SU <- t(apply(S, 1, function(x) proj.tconv(x, V_matrix)))
  } else SU <- S
  if (any(!(apply(SU, 1, is_M_ultrametric, C_matrix=C_matrix)))) stop("BAD PROJECTION")
  # Compute maximal w_min for any x in the intersection of FW(SU) with K.
  w <- max_w_min(K, SU, C_matrix)
  while (is.null(w)) {
    Z <- matrix(rnorm(length(V0), mean=0, sd=sigma_s*w_min_v0), nrow=n)
    S <- V0 + Z
    SU <- t(apply(S, 1, function(x) proj.tconv(x, V_matrix)))
    if (any(!(apply(SU, 1, is_M_ultrametric, C_matrix=C_matrix)))) stop("BAD PROJECTION")
    w <- max_w_min(K, SU, C_matrix)
    sigma_s <- sigma_s*0.97
  }
  w_min <- w_min_C(w, C_matrix)
  return(list(S=SU, w=w, w_min=w_min))
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

# Sample indices of cones in Nested Set fan.
K_x <- sample(seq_len(length(K_NS)), num_samples, replace = TRUE)
# Generate points in those cones.
X <- t(sapply(K_x, function(k) sample_cone(K_NS[[k]])))

S_x <- list()
w_mat <- matrix(0, nrow=num_samples, ncol=length(E_set))
w_min_vec <- c()
for (i in 1:num_samples){
  v0 <- as.numeric(X[i,])
  K <- K_NS[[K_x[i]]]
  Sample.from.U_M <- Generate.Sample.from.U_M(v0, K, sample_size, sigma_s=3)
  S_x[[length(S_x)+1]] <- Sample.from.U_M$S
  w_mat[i,] <- Sample.from.U_M$w
  w_min_vec <- c(w_min_vec, Sample.from.U_M$w_min)
  print(paste0(i,"/", num_samples, " Samples Generated"))
}

save(K_x, X=X, S_x, w_mat, w_min_vec, file = "Data/M(Running)_S1000_N25.rda")
