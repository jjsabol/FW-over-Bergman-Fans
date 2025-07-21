source("Scripts/Functions.R")

# Define M via a list of bases.
EL <- matrix(c(1,2, 2,3, 3,4, 4,1, 1,5, 4,5, 3,6, 4,6), 8, 2, byrow=TRUE)
G <- igraph::graph_from_edgelist(EL, directed=FALSE)
plot(G, edge.label = if (is.null(E(G)$name)) as_ids(E(G)) else E(G)$name)
B_list <- enumerate_sp_trees(G)
# Do you want to use max-plus instead of default (min-plus)?
max_plus <- TRUE

#### Compute all M Information ####

if (TRUE){
  E_set <- sort(unique(unlist(B_list))) # Ground set.
  B_matrix <- do.call(rbind, lapply(B_list, char_vect, E=E_set))
  F_matrix <- proper_flats(B_matrix)
  # Cyclic Bergman Fan, Bergman Fan, and Cycles.
  CBF <- cyclic_Bergman(B_list, silent = FALSE)
  K_Cycl <- lapply(CBF$W, rbind) # Maintains matrix structure in case cones is 1-D
  C_matrix <- do.call(rbind, lapply(CBF$circuits, char_vect, E=E_set))
  F_df <- flat_data(F_matrix, B_matrix, C_matrix)
  # Connected flats for a minimal building set (G_min)
  Rays_NS <- data.frame.to_matrix(F_df[connected==TRUE,1:(length(E_set))])
  G_min <- rbind(rep(1, length(E_set)), Rays_NS)
  # Build the set of maximal cones in the Bergman/Nested Set Fans.
  K_Berg <- list()
  K_Nest <- list()
  i <- 0
  for (partition in CBF$BF_cones){ # For each cone in the Bergman Fan.
    i <- i + 1
    cone_indices <- which(CBF$BF_map == i)
    all_rays <- unique(do.call(rbind, K_Cycl[cone_indices]))
    if (nrow(all_rays) > 1) { 
      BF_cone <- rcdd::redundant(cbind(0,0,all_rays))$output[,-(1:2)]
    } else BF_cone <- all_rays
    K_Berg[[i]] <- BF_cone
    # Compute nested sets and associated cones.
    p_sets <- unlist(strsplit(partition, split = "\\|"))
    p_matrix <- t(sapply(p_sets, char_vect, E=E_set))
    r <- nrow(p_matrix)
    big_pi <- c()
    for (f in 1:nrow(Rays_NS)){
      flat_vec <- as.numeric(Rays_NS[f,])
      s <- max(B_matrix %*% flat_vec) # Rank of the flat.
      combos <- comboGeneral(r, s)
      for (j in 1:nrow(combos)){
        combo <- p_matrix[combos[j,],]
        if (class(combo)[1] != "numeric") combo <- as.numeric(colsums(combo))
        if (all(flat_vec == combo)) { # If flat is a union of s circuits.
          big_pi <- c(big_pi, f) # Add the flat (by its index).
          break
        }
      }
    }
    if (length(big_pi) == r-1) { K_Nest[[length(K_Nest)+1]] <- BF_cone
    } else {
      F_pi <- Rays_NS[big_pi,]
      max_nested <- max_nested_sets(F_pi, G_min)
      for (set in max_nested) {
        K_Nest[[length(K_Nest)+1]] <- F_pi[set,]
      }
    }
  }
  names(K_Berg) <- CBF$BF_cones
  print(paste0("Cyclic Cones: ", length(K_Cycl), "; Nested Set Cones: ", length(K_Nest),
    "; Bergman Cones: ", length(K_Berg)))
}
if (max_plus){
  # Convert cones to Max-Plus
  K_BF <- lapply(K_Berg, function(W) -W)
  K_NS <- lapply(K_Nest, function(W) -W)
  K_CF <- lapply(K_Cycl, function(W) -W)
  V_matrix <- log(1-F_matrix)
} else {
  K_BF <- K_Berg
  K_NS <- K_Nest
  K_CF <- K_Cycl
  V_matrix <- -log(1-F_matrix)
}

#### Mappings between Cones ####

# One-to-Many (will be lists)
BF_CF <- sapply(K_BF, function(w) which_cone(colmeans(w), K_CF))
BF_NS <- sapply(K_BF, function(w) which_cone(colmeans(w), K_NS))
NS_CF <- sapply(K_NS, function(w) which_cone(colmeans(w), K_CF))

# Many-to-One (should be vectors)
NS_BF <- sapply(K_NS, function(w) which_cone(colmeans(w), K_BF))
CF_BF <- sapply(K_CF, function(w) which_cone(colmeans(w), K_BF))
CF_NS <- sapply(K_CF, function(w) which_cone(colmeans(w), K_NS))

# Any not get mapped? (checks)
if(any(sapply(BF_CF, length) == 0)) print("Bad BF_CF Mapping")
if(any(sapply(BF_NS, length) == 0)) print("Bad BF_NS Mapping")
if(any(sapply(NS_BF, length) == 0)) print("Bad NS_BF Mapping")
if(any(sapply(NS_CF, length) == 0)) print("Bad NS_CF Mapping")
if(any(sapply(CF_BF, length) == 0)) print("Bad CF_BF Mapping")
if(any(sapply(CF_NS, length) == 0)) print("Bad CF_NS Mapping")

# Mapping maximal cones to topological codes.
cone_codes <- sapply(K_NS, topological_profile, C_matrix=C_matrix)
if(any(duplicated(cone_codes))) print("DUPLICATED CONE CODES")
# Note: Only performs for relative interior of max cones, not boundary points.

# Pre-Image H-spaces associated to each NESTED SET CONE.
pre_images_K <- lapply(K_NS, pre_K)
H_K <- rep(1:length(K_NS), times = sapply(pre_images_K, length))
H_NS <- purrr::flatten(pre_images_K)

#### Save M Info ####

fn <- paste0("Data/Matroid_Info.rda")
save(E_set, B_matrix, C_matrix, F_matrix, V_matrix, F_df, K_BF, K_NS, K_CF,
  BF_NS, BF_CF, NS_BF, NS_CF, CF_BF, CF_NS, cone_codes, H_K, H_NS, file = fn)

#### Additional Examples for Specifying M ####

# LINEAR MATROIDS
A <- matrix(c(1,0,0, 1,0,0, 0,1,0, 0,1,0, 1,0,-1, 0,1,-1), nrow=3) # Example 2.1 (3x6)
A <- matrix(cbind(diag(4), matrix(c(0,1,0,2,3,1,0,1,1,2,1,0),4,3)), 4, 7) # Example 2.4 (4x7)
A <- matrix(c(1,0,0,0, 1,1,0,0, 1,0,1,0, 1,1,1,0, 1,0,0,1, 1,1,0,1, 1,0,1,1, 1,1,1,1), nrow=4) # Example 4.1 (4x8) Unit Cube
A <- matrix(c(rep(1,16), rep(0:1,each=8), rep(0:1,each=4,times=2), rep(0:1,each=2,times=4), rep(0:1,each=1,times=8)), nrow=5, byrow=TRUE) # Example 4.1 (5x16) Unit Cube
A <- matrix(c(1,0,0,0, 1,0,1,0, 0,1,0,0, 0,1,-1,0, 0,1,0,-1, 0,1,-2,0, 0,1,-1,-1, 0,1,0,-2, 0,1,-3,0, 0,1,-2,-1, 0,1,-1,-2, 0,1,0,-3), 4) # Example 4.2 (4x13)
# Enumerate the bases.
B_list <- reverse_search(A)




# GRAPHICAL MATROIDS
# The following graphs are examples from the original Winter's paper:
# "An algorithm for the enumeration of spanning trees"
G <- igraph::graph_from_data_frame(d = data.frame(from = c(3,2,2,1,1), 
  to = c(4,4,3,3,2), name = seq(5)), directed = FALSE, vertices = seq(4))
G <- igraph::graph_from_data_frame(d = data.frame(from = c(4,2,3,3,2,1), 
  to = c(5,5,5,4,3,2), name = seq(6)), directed = FALSE, vertices = seq(5))
G <- igraph::graph_from_data_frame(d = data.frame(from = c(1,1,1,2,3,4), 
  to = c(2,4,3,3,4,5), name = seq(6)), directed = FALSE, vertices = seq(5))
# More generic graphs.
G <- igraph::make_full_graph(4)
G <- igraph::make_ring(5)
# Random Graphs.
num_nodes <- 5; num_edges <- 7
G <- igraph::sample_gnm(n=num_nodes, m=num_edges)
while(any(igraph::degree(G) == 1) | components(G)$no > 1){ # Ensure connected.
  G <- igraph::sample_gnm(n=num_nodes, m=num_edges)
}
# View the graph.
E(G)$name <- LETTERS[seq_len(ecount(G))]
plot(G, edge.label = if (is.null(E(G)$name)) as_ids(E(G)) else E(G)$name)
# Enumerate the bases.
B_list <- enumerate_sp_trees(G)
# Kirchhoff's Determinant can be used to check the number of bases.
Matrix::det(igraph::laplacian_matrix(G)[-1,-1])
# Every graphical matroid is linear.
A <- get_incidence_matrix(G, directed = TRUE)
if (pracma::Rank(A) < nrow(A)) {
  A <- pracma::rref(A)
  A <- A[which(Rfast::rowsums(A==0) != ncol(A)),]
}
B_list <- reverse_search(A)




# UNIFORM MATROIDS
n <- 4; r <- 2
n_combo <- RcppAlgos::comboGeneral(seq_len(n), r)
B_list <- split(n_combo, seq_len(nrow(n_combo)))
# Or alternatively...
A <- round(10*matrix(runif(n*r), r, n))
while (pracma::Rank(A) < r) A <- round(10*matrix(runif(n*r), r, n))
B_list <- reverse_search(A)





# EXAMPLES FROM FEICHTNER/STURMFELS PAPER
# Example 1.2/3.4 of Feichtner and Sturmfels
K4 <- igraph::make_full_graph(4)
G <- igraph::delete_edges(K4, 1)
E(G)$name <- c(1,2,4,3,5)
plot(G, edge.label = if (is.null(E(G)$name)) as_ids(E(G)) else E(G)$name)
B_list <- enumerate_sp_trees(G)

# Example 2.2 of Feichtner and Sturmfels
B_list <- split(RcppAlgos::comboGeneral(4, 2), seq_len(choose(4,2)))

# Example 2.8 of Feichtner and Sturmfels
U_46 <- RcppAlgos::permuteGeneral(0:1, freqs = c(2,4))
non_bases <- t(sapply(c("1,2,3,4", "1,3,5,6", "2,4,5,6"), char_vect, E=seq(6)))
B_df <- fsetdiff(data.table(U_46), data.table(non_bases))
B_list <- split(t(apply(B_df, 1, function(x) which(x==1))), seq_len(nrow(B_df)))
# We can also represent this matroid as a linear matroid.
L_46 <- RcppAlgos::comboGeneral(6, 4)
A <- round(100*matrix(runif(6*4), 4, 6)) # Start with a uniform matroid.
while (min(apply(L_46, 1, function(x) pracma::Rank(A[,x]))) < 4) A <- round(100*matrix(runif(6*4), 4, 6))
A4 <- A[,1] + A[,2] + A[,3] # Now create dependencies.
A5 <- A[,1] + A[,3] + A[,6]
A[,4] <- A4
A[,5] <- A5
L_46[which(apply(L_46, 1, function(x) pracma::Rank(A[,x])) == 3),] # Check
B_list <- reverse_search(A)

# Example 5.7 of Feichtner and Sturmfels (R_10)
A <- cbind(Rfast::Diag.matrix(5, 1), matrix(c(-1,1,0,0,1, 1,-1,1,0,0, 0,1,-1,1,0, 0,0,1,-1,1, 1,0,0,1,-1), 5, 5))
B_list <- reverse_search(A)

# Example 5.8 of Feichtner and Sturmfels (Cographic M(K5))
K5 <- igraph::make_full_graph(5)
sp_list <- enumerate_sp_trees(K5)
B_list <- lapply(sp_list, function(x) setdiff(seq_len(10), x))

# Example 5.9 of Feichtner and Sturmfels (3D Unit Cube)
faces <- matrix(c(1,2,3,4, 5,6,7,8, 1,2,5,6, 3,4,7,8, 1,3,5,7, 2,4,6,8), 6, 4, byrow=TRUE)
nonfaces <- fsetdiff(data.table(comboGeneral(8, 4)), data.table(faces))
B_list <- split(data.frame.to_matrix(nonfaces), seq_len(nrow(nonfaces))) 

# Example 5.9 (4D Unit Cube)
vert <- permuteGeneral(0:1, 4, repetition = TRUE)
A <- rbind(t(vert), 1)
#B_list <- reverse_search(A) # SLOW
combos <- comboGeneral(16, 5)
combo_ind <- apply(combos, 1, function(x) Matrix::rankMatrix(A[,x])[1] == 5)
B_list <- split(combos[which(combo_ind),], seq_len(sum(combo_ind)))
