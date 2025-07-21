# Required packages.
library(data.table)
library(Rfast)
library(RcppAlgos)
library(igraph)
library(fpCompare)
library(rcdd)
library(lpSolve)
library(ggplot2)
library(Rglpk)
library(purrr)

#### General Functions ####

# Returns the characteristic vector based on ground set E.
char_vect <- function(E,v){
  # Characteristic vector given ground set E.
  if (class(v) == "character") v <- unlist(strsplit(v,","))
  if (class(E) == "numeric") v <- as.numeric(v)
  return((E %in% v)+0)
}
# Hamming distance between two vectors.
hamming <- function(x,y) {
  # Get Hamming distance from x to (all rows of) y.
  if ("data.frame" %in% class(y) | "matrix" %in% class(y)) {
    return(Rfast::rowsums(Rfast::eachrow(y, x, oper = "-") %!=% 0))
  } else return(sum(x != y))
}
# Incidence matrix of a graph.
get_incidence_matrix <- function(g, directed = FALSE) {
  # Get incidence matrix from a graph.
  n <- vcount(g)
  m <- ecount(g)
  M <- matrix(0, nrow = n, ncol = m)
  edges <- as_edgelist(g, names = FALSE)
  for (j in seq_len(m)) {
    v1 <- edges[j, 1]
    v2 <- edges[j, 2]
    M[v1, j] <- 1
    M[v2, j] <- 1 - 2*directed
  }
  rownames(M) <- V(g)$name
  colnames(M) <- paste0("e", seq_len(m))
  return(M)
}
# Tropical projection onto convex hull of V.
proj.tconv <- function(x, V, max.plus=TRUE){
  # Project onto the tropical convex hull.
  # Description: Equivalent to TML::project_pi(), but without for-loops.
  # Inputs: Point to project (x0), onto the polytope defined by vertices in V.
  # Output: The projected point.
  # Dependencies: None
  # Ex: W <- matrix(c(0,0,0, 0,2,3, 0,5,1), 3, 3); v <- c(0,7,-1)
  if (class(x)[1] != "numeric") x <- as.numeric(x)
  if (class(V)[1] != "matrix") V <- data.frame.to_matrix(V)
  if (max.plus){
    lambda <- rowMins(eachrow(-V, x, oper = "+"), value=T)
    pi.V.x <- rowMaxs(eachrow(Rfast::transpose(V), lambda, oper = "+"), value=T)
  } else{
    lambda <- rowMaxs(eachrow(-V, x, oper = "+"), value=T)
    pi.V.x <- rowMins(eachrow(Rfast::transpose(V), lambda, oper = "+"), value=T)
  }
  return(pi.V.x)
}
# Various ways to normalize a vector/matrix.
norm.i <- function(X, i=0, center=FALSE, diag=FALSE){
  # Dependencies: Rfast, data.table
  # Notes: If i=0, will provide the "canonical coordinates" in the sense of (2), Section 5.2,
  # otherwise, specifying i as an integer will set the coordinates of that dimension to zero.
  # Notes: If center==TRUE, then each row sum will be zero as in (1).
  if (class(X)[1] == "numeric"){
    if (center==TRUE){ X.norm <- X - mean(X)
    } else if (i==0){ X.norm <- X - min(X)
    } else X.norm <- X - X[i]
  } else{
    cn <- colnames(X)
    if (class(X)[1] != "matrix") { X.mat <- Rfast::transpose(Rfast::data.frame.to_matrix(X))
    } else {
      X.mat <- try(Rfast::transpose(X), silent = TRUE)
      if ("try-error" %in% class(X.mat)) X.mat <- try(Rfast::transpose(matrix(unlist(X), nrow(X), ncol(X))), silent = TRUE)
    }
    if (center==TRUE){ X.norm <- Rfast::transpose(Rfast::eachrow(X.mat, colMeans(X.mat, na.rm = TRUE), oper = "-"))
    } else if (diag==TRUE & nrow(X.mat) == ncol(X.mat)){ X.norm <- X.mat - rep_row(diag(X.mat), 3)
    } else if (i==0){ X.norm <- Rfast::transpose(Rfast::eachrow(X.mat, Rfast::colMins(X.mat, value=TRUE), oper = "-"))
    } else X.norm <- Rfast::transpose(Rfast::eachrow(X.mat, X.mat[i,], oper = "-"))
    colnames(X.norm) <- cn
  }
  if ("data.table" %in% class(X)) X.norm <- data.table::as.data.table(X.norm)
  return(X.norm)
}

# Cone Inclusion
# Sample from the relative interior of a cone.
sample_cone <- function(K, ray_weights=NULL, unit_len=TRUE){
  r <- nrow(K)
  if (is.null(ray_weights)) { alpha <- rep(1, r) # Uniform weighting.
  } else if (length(ray_weights) == r) { alpha <- ray_weights
  } else stop("Invalid weights.")
  gamma_samples <- rgamma(r, shape = alpha, rate = 1)
  lambda <- gamma_samples / sum(gamma_samples)
  s <- rexp(1, rate = 1)
  x <- as.numeric(s * (t(K) %*% as.vector(lambda)))
  if (unit_len) x <- x/max(abs(x))
  return(x)
}
# Test (via LP) if a vector is in a particular cone.
in_cone <- function(x, W, add1=TRUE){
  if (add1) W <- rbind(1, W)
  lp_test <- lp("min", objective.in = rep(0, nrow(W)), const.mat = t(W), 
    const.dir = rep("==", ncol(W)), const.rhs = x)
  return(lp_test$status == 0) 
}
# Which cones (among a set of cones) is the vector in?
which_cone <- function(x, W, add1=TRUE){
  x <- x - min(x) # Ensure non-negative x.
  if (diff(min_max(x)) %>>% 0){
    x <- x/min(x[which(x %>>% 0)]) # This helps avoid numerical issues with lpSolve.
    return(which(sapply(W, in_cone, x=x, add1=add1)))
  } else return(seq_len(length(W)))
}
# An alternative method for determining cone inclusion. Maps between cone codes.
recover_K <- function(x){
  tx <- topological_profile(x, C_matrix)
  j <- match(tx, cone_codes)
  if (is.na(j)){
    if (tx %in% names(boundary_cases)){ # Check if it's in already computed boundary instances.
      I <- boundary_cases[[tx]]$Bergman
      J <- boundary_cases[[tx]]$Nested
      K <- boundary_cases[[tx]]$Cyclic
    } else { # Otherwise, add it to the boundary list.
      K <- which_cone(x, K_CF) # Test all cone membership in the finest structure.
      J <- CF_NS[K] # Translate them into membership of the coarser structures.
      I <- CF_BF[K]
      boundary_cases[[tx]] <<- list(Bergman = I, Nested = J, Cyclic = K)
    }
  } else {
    J <- j
    I <- NS_BF[j]
    possible_K <- unlist(NS_CF[j])
    if (length(possible_K) > 1) { K <- possible_K[which_cone(x, K_CF[possible_K])]
    } else K <- possible_K
  }
  BF_x <- char_vect(seq_len(max(NS_BF)), I)
  NS_x <- char_vect(seq_len(length(NS_CF)), J)
  CF_x <- char_vect(seq_len(length(CF_NS)), K)
  BF_x <- BF_x/sum(BF_x)
  NS_x <- NS_x/sum(NS_x)
  CF_x <- CF_x/sum(CF_x)
  #return(data.table(Bergman=i, Nested=j, Cyclic=k))
  return(list(Bergman=BF_x, Nested=NS_x, Cyclic=CF_x))
}
# Similar to "in_cone" but using H-representation to solve LP.
in_Hspace <- function(x, K) all(K %*% x %>=% 0)
# Which polyhedron (among a set of them) is the vector in?
which_Hspace <- function(x, K_set) which(sapply(K_set, function(K) in_Hspace(x,K)))
# Compute the pre-image (under tropical projection) of maximal cone K
pre_K <- function(K){
  E <- seq(ncol(K))
  M <- list()
  rep_j <- c(1)
  for (j in E[-1]){
    if (!any(apply(K[,1:(j-1),drop=FALSE], 2, function(x) all(x == K[,j])))) rep_j <- c(rep_j, j)
  }
  for (j in rep_j){
    j_lt_k <- c()
    k_gt_j <- c()
    K_j <- K[,j]
    for (k in E[-j]){
      K_k <- K[,k]
      if (all(K_j %<=% K_k)) {
        if (any(K_j %<<% K_k)) {
          j_lt_k <- c(j_lt_k, j)
          k_gt_j <- c(k_gt_j, k)
        } 
        if (all(K_j %==% K_k)){
          j_lt_k <- c(j_lt_k, k)
        }
      }
    }
    j_lt_k <- unique(j_lt_k)
    k_gt_j <- unique(k_gt_j)
    if (length(j_lt_k)>0 & length(k_gt_j)>0){
      temp_M <- list()
      for (j in j_lt_k){
        Mi <- matrix(0, length(k_gt_j), length(E))
        r <- 0
        for (k in k_gt_j){
          r <- r+1
          Mi[r,] <- char_vect(E, k) - char_vect(E, j)
        }
        if (length(M) > 0){
          new_M <- M
          for (i in 1:length(new_M)){
            new_M[[i]] <- rbind(new_M[[i]], Mi)
            temp_M[[length(temp_M)+1]] <- rbind(new_M[[i]], Mi)
          }
        } else temp_M[[length(temp_M)+1]] <- Mi
      }
      M <- temp_M
    }
  }
  M_out <- lapply(M, unique)
  return(M_out)
}

#### Algorithms ####

# Algorithmic implementations to get bases.
# Reverse search of Avis/Fukuda.
reverse_search <- function(A){
  # Basis enumeration of columns of matrix A.
  # Subroutines:
  is_basis <- function(A, cols){
    A_rank <- Matrix::qr(A)$rank
    submatrix <- A[, cols, drop = FALSE]
    return(Matrix::qr(submatrix)$rank == A_rank)
  }
  adj_bases <- function(A, B){
    # Requires is_basis
    n <- ncol(A)
    if (class(B) == "numeric") B <- which(B == 1)
    r <- length(B)
    unused <- setdiff(1:n, B)
    adj_list <- list()
    for (i in seq_along(B)) {
      for (e in unused) {
        new_cols <- sort(c(B[-i], e))
        if (is_basis(A, new_cols)) adj_list[[length(adj_list)+1]] <- new_cols
      }
    }
    return(adj_list)
  }
  parent_B <- function(A, B){
    # Requires adj_bases
    lex_order <- function(L){
      max_len <- max(sapply(L, length))
      padded <- t(sapply(L, function(x) c(x, rep(Inf, max_len - length(x)))))
      order_L <- do.call(order, as.data.frame(padded))
      return(order_L)
    }
    nbrs <- adj_bases(A, B)
    if (length(nbrs) == 0) return(NULL)
    nbrs[[length(nbrs)+1]] <- B
    lex_ord_nbrs <- lex_order(nbrs)
    parent_of_B <- nbrs[[lex_ord_nbrs[1]]]
    if (sum(parent_of_B != B) > 0) { return(parent_of_B)
    } else return(NULL)
  }
  # Main:
  n <- ncol(A)
  r <- Matrix::qr(A)$rank
  visited <- new.env(hash = TRUE)
  bases <- list()
  # Start with lex smallest basis
  possible_bases <- comboGeneral(n, r)
  for (i in 1:nrow(possible_bases)) {
    if (is_basis(A, possible_bases[i, ])) {
      root <- possible_bases[i, ]
      break
    }
  }
  dfs <- function(basis) {
    key <- paste0(basis, collapse = ",")
    if (!exists(key, envir = visited)) {
      visited[[key]] <<- TRUE
      bases[[length(bases) + 1]] <<- basis
      for (nbr in adj_bases(A, basis)) {
        if (identical(parent_B(A, nbr), basis)) {
          #print(paste0(length(visited), "---", paste0(nbr, collapse = ",")))
          dfs(nbr)
        }
      }
    }
  }
  dfs(root)
  return(bases)
}
# Spanning trees via Winter.
enumerate_sp_trees <- function(G){
  # Basis enumeration of spanning trees in graph using Winter's algorithm.
  if (!is_connected(G)) stop("G needs to be a connected graph.")
  # SUBROUTINES
  proc_label <- function(G){
    v_labels <- rep(0, vcount(G))
    i <- 1
    v_labels[which(igraph::degree(G) == igraph::max_degree(G))[1]] <- i
    if (vcount(G) > 1) {
      for (each_node in 2:vcount(G)){
        i <- i + 1
        labelled <- which(v_labels != 0)
        adj_vert <- unique(unlist(neighborhood(G, nodes = labelled)))
        v_labels[adj_vert[which(v_labels[adj_vert] == 0)[1]]] <- i
      }
    }
    return(v_labels)
  }
  proc_initialize <- function(G){
    E <- list()
    for (j in 2:n){
      E[[j]] <- list()
      for (i in igraph::neighbors(G, j)) E[[j]][[i]] <- igraph::E(G)[j %--% i]$name
    }
    EE <- list()
    for (j in 2:n) EE[[j]] <- which(sapply(E[[j]], length) > 0)
    return(list(E=E, EE=EE))
  }
  MM <- function(j) {
    EE_j <- EE[[j]]
    if (length(EE_j) > 0) EE_j <- EE_j[EE_j < j]
    if (length(EE_j) == 0) return(0)
    return(max(EE_j))
  }
  F_ji <- function(j,i) {
    E_ji <- E[[j]][[i]]
    if (length(E_ji) > 0) { return(head(E_ji,1))
    } else return(0)
  }
  proc_output <- function(){
    cp <- E[[n]][[r[n]]]
    for (i in 1:(n-2)){
      if ("data.frame" %in% class(cp)) {
        cp[, i := seq_len(nrow(cp))]
        cp2 <- data.table::CJ(i = cp$i, E[[n-i]][[r[n-i]]])
        cp <- merge.data.table(cp, cp2, by = "i", all.y = TRUE)
        colnames(cp) <- c("i", paste0("k.", seq(ncol(cp)-1)))
      } else cp <- data.table::CJ(n=cp, E[[n-i]][[r[n-i]]])
    }
    cp[, i := NULL]
    sp_trees <<- rbindlist(list(sp_trees, cp))
    return(cp)
  }
  proc_contract <- function(){
    nk <- n-k # Node to be removed.
    rk <- MM(nk) # Node to be merged.
    while (rk != 0){
      r[nk] <<- rk
      if ((n-k) == 2) {
        z <- length(E[[n]][[r[n]]])
        for (i in (n-1):2) z <- min(z, length(E[[n-1]][[r[n-1]]]))
        if (z > 0) { return(proc_output())
        } else return()
      } else {
        EE_nk <- EE[[nk]] # Set of nodes adjacent to removed node.
        EE_set <- EE_nk[EE_nk < rk]
        for (i in rev(EE_set)){
          if (length(E[[rk]][[i]]) > 0){
            E_ji <- c(E[[rk]][[i]], E[[nk]][[i]][E[[nk]][[i]] %notin% E[[rk]][[i]]])
            E[[rk]][[i]] <<- E_ji
          } else {
            E[[rk]][[i]] <<- E[[nk]][[i]]
            if (i %notin% EE[[rk]]) EE[[rk]] <<- c(i, EE[[rk]])
          }
        }
        k <<- k + 1
        proc_contract()
        k <<- k - 1
        rpp <- 0
        EE_nk <- EE[[nk]]
        EE_set <- EE_nk[EE_nk < rk]
        for (i in EE_set){
          E[[rk]][[i]] <<- E[[rk]][[i]][E[[rk]][[i]] %notin% E[[nk]][[i]]]
          if (F_ji(nk,i) == F_ji(rk,i)){
            EE[[rk]] <<- EE[[rk]][EE[[rk]] != i]
          }
          if (i > rpp) rpp <- i
        }
        rk <- rpp
      }
    }
  }
  # BEGIN WINTERS
  n <- igraph::vcount(G)
  old_v_labs <- if(is.null(V(G)$name)) as_ids(V(G)) else V(G)$name
  if (n == 1) return(old_v_labs)
  old_e_labs <- if(is.null(E(G)$name)) as_ids(E(G)) else E(G)$name
  V(G)$name <- proc_label(G)
  E(G)$name <- old_e_labs
  init <- proc_initialize(G)
  E <- init$E
  EE <- init$EE
  k <- 0
  r <- rep(0, n)
  sp_trees <- data.table()
  proc_contract()
  sp_list <- split(as.matrix(sp_trees), seq_len(nrow(sp_trees)))
  #sp_list <- lapply(sp_list, sort)
  return(sp_list)
}

# Compute fan structures.
# Cyclic Bergman via Rincon.
cyclic_Bergman <- function(B_list, silent=FALSE){
  E <- sort(unique(unlist(B_list)))
  n <- length(E)
  # Subroutines.
  F_k <- function(k){
    # Requires global access to B and B_list.
    I_k <- c()
    for (i in B){
      Bik <- c(setdiff(B,i), k)
      ind <- any(sapply(B_list, function(basis) all(Bik %in% basis)))
      if (ind) I_k <- c(I_k,i)
    }
    return(sort(I_k))
  }
  L_smallest <- function(x,L) x[which.min(match(x, L))[1]]
  extend_L <- function(L,b) {
    r <- length(L)
    if (r == 0) return(b)
    result <- vector("list", r + 1)  # We will get r+1 permutations
    for (i in 0:r) result[[i + 1]] <- append(L, b, after = i)
    return(result)
  }
  Q_cones <- function(p,L,F_nB){
    # The following could theoretically pull from globals (except B for some reason).
    n <- length(p)
    B <- which(is.na(p)) # This specific line seems to be required for cones to be calculated correctly.
    nB <- which(!is.na(p))
    #for (k in nB) F_nB[[k]] <- F_k(k)
    r <- length(B)
    # BEGIN
    Qb <- matrix(0, r, n)
    Wb <- matrix(0, r, n)
    for (i in 1:r) Qb[i, c(B[i], which(p == B[i]))] <- 1
    # OPTIONAL. Associates the cyclic cones to the Bergman cones.
    partitions <- sort(apply(Qb, 1, function(x) paste0(E[which(x == 1)], collapse = ",")))
    partition_string <- paste0(partitions, collapse = "|")
    associated_cone <- which(BF_cones == partition_string)
    if (length(associated_cone) == 0) { 
      BF_cones <<- c(BF_cones, partition_string)
      associated_cone <- length(BF_cones)
    }
    BF_map <<- c(BF_map, associated_cone)
    # END OPTIONAL.
    edgelist <- lapply(seq_len(length(L)-1), function(i) L[i:(i+1)])
    for (cc in B[B %notin% p]){
      b <- c()
      for (k in nB) if (any(cc == F_nB[[k]])) b <- c(b, p[k])
      if (length(b) > 1) b <- b[which.max(match(b, L))]
      edgelist[[length(edgelist)+1]] <- c(b, cc)
    }
    edgelist <- matrix(unlist(edgelist), ncol = 2, byrow = TRUE)
    Tb <- igraph::graph_from_edgelist(el = edgelist)
    Tb <- delete_vertices(Tb, V(Tb)[igraph::degree(Tb)==0])
    for (i in 1:length(B)){
      reachable <- igraph::subcomponent(Tb, v=i, mode="out")
      if(length(reachable) > 1){ Wb[i,] <- Rfast::colsums(Qb[reachable,])
      } else Wb[i,] <- Qb[reachable,]
    }
    Wb <- Wb[which(Rfast::rowsums(Wb) < n),] # Mods out (1,1...,1).
    return(Wb)
  }
  proc_Pref <- function(k,p,L, verbose=FALSE){
    # Requires: L_smallest, Q, extend_L
    # Needs global access to F_nB, nB
    if (k == "end"){
      counter <<- counter + 1
      if (verbose) print(paste0(counter, "-  p: (", paste0(p, collapse = ","), "); L: (", paste0(L, collapse = ","), ")"))
      W_list[[counter]] <<- Q_cones(p,L,F_nB)
      p_list <<- c(p_list, paste0(p, collapse = ","))
      L_list <<- c(L_list, paste0(L, collapse = ","))
      return()
    } else{
      Fk <- F_nB[[k]]
      nB_gt_k <- nB[nB > k]
      nB_lt_k <- nB[nB < k]
      L_int_Fk <- intersect(L, Fk)
      if (length(L_int_Fk) > 0){
        p[k] <- if(length(L_int_Fk) == 1) L_int_Fk else L_smallest(L_int_Fk, L)
        kp <- if(length(nB_gt_k) == 0) "end" else nB_gt_k[1]
        proc_Pref(kp, p, L, verbose)
      }
      for (b in setdiff(Fk[Fk < k], L)){
        for (Lp in extend_L(L,b)){
          if (b == L_smallest(Fk, Lp)){ #<- RINCON DID NOT INCLUDE THIS CHECK IN HIS PSEUDOCODE!
            inadmissible <- FALSE
            for (l in nB_lt_k){
              inadmissible <- any(b == F_nB[[l]]) & match(b, Lp) < match(p[l], Lp)
              if (inadmissible) break
            }
            if (!inadmissible){
              p[k] <- b
              kp <- if(length(nB_gt_k) == 0) "end" else nB_gt_k[1]
              proc_Pref(kp, p, Lp, verbose)
            }
          }
        }
      }
    }
  }
  # Outer Initialization
  counter <- 0
  W_list <- list()
  p_list <- c()
  L_list <- c()
  C_list <- c()
  BF_cones <- c()
  BF_map <- c()
  F_nB <- list()
  for (B in B_list) {
    nB <- setdiff(E, B)
    for (k in nB) {
      Fk <- F_k(k)
      C_str <- paste0(sort(c(Fk,k)), collapse = ",")
      if (length(C_list) == 0) { C_list <- C_str
      } else if (!any(C_str == C_list)) C_list <- c(C_list, C_str)
      if (class(B) == "character"){
        k <- match(k, E)
        Fk <- match(Fk, E)
      }
      F_nB[[k]] <- Fk
    }
    if (class(B) == "character"){
      B <- match(B, E)
      nB <- match(nB, E)
    }
    # Inner Initialization
    p <- rep(NA, n)
    L <- c()
    k <- nB[1]
    proc_Pref(k, p, L, verbose = !silent)
  }
  return(list(E=E, W=W_list, p=p_list, L=L_list, BF_cones=BF_cones, BF_map=BF_map, circuits=C_list))
}
# Nested Sets via Feichtner/Sturmfels.
max_nested_sets <- function(W_pi, G){
  incomparable <- function(wi, wj) !(any(min_max(wi - wj) == 0))
  in_build_set <- function(w, G) any(apply(Rfast::eachrow(G, w, oper = "-"), 1, function(x) all(min_max(x) == 0)))
  N <- seq_len(nrow(W_pi))
  nested_sets <- split(N, N)
  for (k in 2:nrow(W_pi)){
    bigger_sets <- list()
    for (set in nested_sets){
      for (j in setdiff(N, set)){
        incomp <- sapply(set, function(ii) incomparable(W_pi[ii,], W_pi[j,]))
        w_ij <- sapply(set, function(k) as.numeric(Rfast::colsums(W_pi[c(k, j),]) > 0))
        inBset <- apply(w_ij, 2, function(wij) in_build_set(wij, G))
        if (all((incomp & !inBset) | !incomp)) bigger_sets[[length(bigger_sets)+1]] <- sort(c(set, j))
      }
    }
    bigger_sets <- unique(bigger_sets)
    if (length(bigger_sets) == 0) break
    nested_sets <- unique(bigger_sets)
  }
  return(nested_sets)
}

#### Matroid Functions ####

fundamental_circuit <- function(k,B,B_sets){
  # Fundamental Circuit
  # B_sets can be given as either list or (NOT 0/1) matrix.
  I_k <- k
  for (i in B){
    Bik <- c(setdiff(B,i), k)
    if ("matrix" %in% class(B_sets)){
      ind <- any(apply(B_sets, 1, function(basis) all(Bik %in% basis)))
    } else {
      ind <- any(sapply(B_sets, function(basis) all(Bik %in% basis))) 
    }
    if (ind) I_k <- c(I_k,i)
  }
  return(sort(I_k))
}
circuits <- function(B_matrix){
  # Get all circuits (cyclic Bergman does this automatically).
  n <- ncol(B_matrix)
  E <- Rfast::Diag.matrix(n, v=1)
  r <- sum(B_matrix[1,])
  C_matrix <- matrix(1, 1, n)
  for (m in 2:(r+1)){
    c_m <- comboGeneral(n, m=m)
    Cm <- t(apply(c_m, 1, function(x) Rfast::colsums(E[x,])))
    ind <- apply(Cm, 1, function(x) any(c(B_matrix %*% x) == m))
    if (any(!ind)) {
      dep_Cm <- Cm[which(!ind),]
      if ("matrix" %in% class(dep_Cm)) {
        sub_C <- apply(dep_Cm, 1, function(x) any(c(C_matrix %*% x) == rowsums(C_matrix)))
        if (any(!sub_C)) C_matrix <- rbind(C_matrix, dep_Cm[which(!sub_C),])
      } else {
        if (nrow(C_matrix)==1) { C_matrix <- rbind(C_matrix, dep_Cm)
        } else {
          sub_C <- c(C_matrix %*% dep_Cm)[-1] == rowsums(C_matrix)[-1]
          if (!sub_C) C_matrix <- rbind(C_matrix, dep_Cm)
        }
      }
    }
  }
  return(C_matrix[-1,])
}
M_connected <- function(C_matrix){
  # Equivalence relations on the ground set [n], c(M).
  n <- ncol(C_matrix)
  IJ <- comboGeneral(n,2)
  E <- Rfast::Diag.matrix(n, v=1)
  A <- Rfast::Diag.matrix(n, v=0)
  for (k in 1:nrow(IJ)){
    i <- IJ[k,1]
    j <- IJ[k,2]
    e_ij <- colsums(E[IJ[k,],])
    if (max(C_matrix %*% e_ij) == 2) A[i,j] <- 1
  }
  A <- A + t(A)
  G_A <- igraph::graph_from_adjacency_matrix(A)
  igraph::components(G_A)$no == 1
}
# Flat-specific functions.
is_flat <- function(flat_vec, B_matrix, C_matrix=NULL){
  # Is a 0/1 vector a flat given the list of bases or circuits?
  if (!is.null(C_matrix)){
    C_f <- Rfast::eachrow(C_matrix, flat_vec, oper = "-")
    return(sum(Rfast::rowsums(C_f %==% 1) %==% 1) %==% 0)
  } else{
    flat_vec <- as.numeric(flat_vec)
    rF <- max(B_matrix %*% flat_vec)
    for (e in which(flat_vec == 0)){
      F_e <- flat_vec
      F_e[e] <- 1
      rF_e <- max(B_matrix %*% F_e)
      if (rF_e == rF) return (FALSE)
    }
    return (TRUE)
  }
}
proper_flats <- function(B_matrix){
  # Enumerate all proper flats given the set of all circuits.
  n <- ncol(B_matrix)
  F_sets <- data.table()
  for (i in 1:(n-1)){
    n_choose_i <- data.table(permuteGeneral(1:0, freqs = c(i, n-i)))
    F_sets <- rbindlist(list(F_sets, n_choose_i))
  }
  F_mat <- data.frame.to_matrix(F_sets)
  apply(F_mat, 1, function(x) is_flat(x, B_matrix))
  F_sets[, is.flat := apply(data.frame.to_matrix(F_sets), 1, is_flat, B_matrix=B_matrix)]
  return(data.frame.to_matrix(subset(F_sets, is.flat)[,1:n]))
}
flat_connected <- function(flat_vec, B_matrix, upper_interval=FALSE){
  # Given a flat, is it connected?
  flat_vec <- as.numeric(flat_vec)
  B_int_f <- Rfast::eachrow(B_matrix, flat_vec)
  B_int_f_rank <- Rfast::rowsums(B_int_f)
  s <- max(B_int_f_rank)
  if (upper_interval){
    #B_min_f <- (Rfast::eachrow(B_matrix, flat_vec, oper = "-") > 0)+0
    #B_F <- B_int_f[B_int_f_rank == s,][1,]
    #B_BF <- Rfast::eachrow(B_matrix, B_F, oper = "-")
    #M_F <- B_min_f[which(Rfast::rowMins(B_BF, value=TRUE) >= 0),]
    B_min_f <- (Rfast::eachrow(B_matrix, flat_vec, oper = "-") > 0)+0
    M_F <- B_min_f[which(B_int_f_rank == s),] # The CONTRACTION of M to F.
  } else M_F <- B_int_f[which(B_int_f_rank == s),] # The RESTRICTION of M to F.
  if ("matrix" %in% class(M_F)) {
    M_F <- unique(M_F)
    M_F <- M_F[,colsums(M_F) > 0]
    if (class(M_F)[1] == "numeric"){ conn <- sum(M_F) == 1
    } else {
      M_F_circuits <- circuits(M_F)
      if (class(M_F_circuits)[1] == "numeric") { conn <- sum(M_F_circuits) == ncol(M_F)
      } else conn <- M_connected(M_F_circuits)
    }
  } else conn <- sum(M_F) == 1
  return(conn)
}
flat_cyclic <- function(flat_vec, C_matrix){
  # Given a flat, is it cyclic?
  f <- as.numeric(flat_vec)
  C_f <- Rfast::eachrow(C_matrix, f, "-")
  p <- which(rowMaxs(C_f, value=TRUE) <= 0)
  if (length(p) == 0) return(FALSE)
  if (length(p) == 1) return(all(f == as.numeric(C_matrix[p,])))
  return(all(f == (Rfast::colsums(C_matrix[p,]) > 0)+0))
}
flat_data <- function(F_matrix, B_matrix, C_matrix){
  # Characterizes all flats in a data table.
  Flats <- data.table(F_matrix)
  is_singleton <- Rfast::rowsums(F_matrix) == 1
  is_cyclicFlt <- apply(F_matrix, 1, flat_cyclic, C_matrix=C_matrix)
  is_connected <- apply(F_matrix, 1, flat_connected, B_matrix=B_matrix)
  is_upperConn <- apply(F_matrix, 1, flat_connected, B_matrix=B_matrix, upper_interval=TRUE)
  Flats[, is_single := is_singleton]
  Flats[, is_cyclic := is_cyclicFlt]
  Flats[, connected := is_connected]
  Flats[, is_flacet := is_connected & is_upperConn]
  print(paste0("Flats: ", nrow(F_matrix), "; Cyclic/Single: ", sum(is_singleton | is_cyclicFlt), 
    "; Connected: ", sum(is_connected), "; Flacets: ", sum(is_connected & is_upperConn)))
  return(Flats)
}
# M-Ultrametric
M_subDomUltra <- function(w, B_matrix){
  # Project onto the tropical convex hull.
  # Slower than proj.tconv, but you only need bases.
  E <- seq_len(length(w))
  E_matrix <- Rfast::Diag.matrix(length(w), 1)
  Bw <- Rfast::rowsums(Rfast::eachrow(B_matrix, w))
  Bw_min <- which(B_matrix[Rfast::min_max(Bw, index = TRUE)[1],] %!=% 0)
  #Bw_max <- which(B_matrix[Rfast::min_max(Bw, index = TRUE)[2],] %!=% 0)
  B_sets <- t(apply(B_matrix, 1, function(x) which(x == 1)))
  w_M <- w
  for (e in setdiff(E, Bw_min)){
    C0 <- fundamental_circuit(e, Bw_min, B_sets)
    C0e <- setdiff(C0, e)
    w_M[e] <- max(w[C0e])
  }
  return(w_M)
}
is_M_ultrametric <- function(x, C_matrix, max.plus = TRUE){
  # Is point an M-ultrameteric?
  if (max.plus) {
    Cx <- Rfast::eachrow(log(C_matrix), x, oper = "+")
    C_maxs <- Rfast::rowMaxs(Cx, value = TRUE)
    achieve_ext <- Rfast::colsums(Rfast::eachrow(Rfast::transpose(Cx), C_maxs, oper = "-") %>=% 0)
  } else {
    Cx <- Rfast::eachrow(-log(C_matrix), x, oper = "+")
    C_mins <- Rfast::rowMins(Cx, value = TRUE)
    achieve_ext <- Rfast::colsums(Rfast::eachrow(Rfast::transpose(Cx), C_mins, oper = "-") %<=% 0)
  }
  return (min(achieve_ext) >= 2)
}
w_min_C <- function(w, C_matrix){
  # Per the paper.
  Cw <- Rfast::eachrow(log(C_matrix), w, oper = "+")
  Cw_max <- Rfast::rowMaxs(Cw, value = TRUE)
  Tw <- Rfast::eachrow(Rfast::transpose(Cw), Cw_max, oper = "-") %>=% 0
  Cw2 <- Cw + log(1-Rfast::transpose(Tw))
  Cw_2max <- Rfast::rowMaxs(Cw2, value = TRUE)
  w_min <- min(Cw_max - Cw_2max)
  return(w_min)
}
topological_profile <- function(w, C_matrix, max_plus=TRUE, as_code=TRUE){
  if ("matrix" %in% class(w)) w <- Rfast::colmeans(w)
  # Records the argmax of each circuit.
  if (max_plus) {
    Cw <- Rfast::eachrow(log(C_matrix), w, oper = "+")
    Cw_maxs <- Rfast::rowMaxs(Cw, value=TRUE)
    Cw_max2 <- t(Rfast::eachrow(t(Cw), Cw_maxs, oper = "-"))
    argmax <- (Cw_max2 %==% 0)+0
    profile_code <- paste0(apply(argmax, 1, function(x) paste0(which(x==1), collapse = ",")), collapse = "|")
    if (as_code) { return(profile_code)
    } else return(argmax)
  } else {
    Cw <- Rfast::eachrow(-log(C_matrix), w, oper = "+")
    Cw_mins <- Rfast::rowMins(Cw, value=TRUE)
    Cw_min2 <- t(Rfast::eachrow(t(Cw), Cw_mins, oper = "-"))
    argmin <- (Cw_min2 %==% 0)+0
    profile_code <- paste0(apply(argmin, 1, function(x) paste0(which(x==1), collapse = ",")), collapse = "|")
    if (as_code) { return(profile_code)
    } else return(argmin)
  }
}

# Graphs
G_BasisExchange <- function(B_list){
  # Basis Exchange Graph
  r <- length(B_list[[1]])
  d <- rep(NA, choose(length(B_list),2))
  k <- 1
  for (i in 1:(length(B_list)-1)){
    for (j in (i+1):length(B_list)){
      d[k] <- length(intersect(B_list[[i]], B_list[[j]]))
      k <- k+1
    }
  }
  d <- (d==r-1)+0
  BE_mat <- matrix(0, length(B_list), length(B_list))
  BE_mat <- Rfast::lower_tri.assign(BE_mat, d)
  BE_mat <- BE_mat + Rfast::transpose(BE_mat)
  g <- igraph::graph_from_adjacency_matrix(BE_mat, mode = "undirected")
  if (vcount(g) <= 10){
    v_labs <- unlist(lapply(B_list, paste0, collapse = ","))
    V(g)$label <- sapply(v_labs, function(x) paste0("{", x, "}"))
    plot(g, vertex.size=3, vertex.label.cex=0.7, vertex.label.dist=1, vertex.label.degree=pi/2)
  } else plot(g, vertex.size=3, vertex.label=NA)
  return(g)
}
G_LatticeFlats <- function(F_matrix, B_matrix){
  # Lattice of Flats
  flat_ranks <- apply(F_matrix, 1, function(x) max(B_matrix %*% as.numeric(x)))
  flat_ranks <- c(flat_ranks, 0, max(flat_ranks)+1)
  flats <- split(F_matrix, seq_len(nrow(F_matrix)))
  flats <- lapply(flats, function(x) which(x==1))
  flats[[length(flats)+1]] <- c(numeric(0))
  flats[[length(flats)+1]] <- seq_len(ncol(F_matrix))
  # Function to test if x is subset of y
  is_sub_flat <- function(x, y) all(x %in% y)
  # Find cover relations
  edges <- list()
  for (i in 1:(length(flats)-1)) {
    for (j in (i+1):length(flats)) {
      rank_i <- flat_ranks[i]
      rank_j <- flat_ranks[j]
      flat_i <- flats[[i]]
      flat_j <- flats[[j]]
      if (rank_i - rank_j == 1){
        if (rank_j == 0) { edges[[length(edges)+1]] <- c(j,i)
        } else if (is_sub_flat(flat_j, flat_i)) edges[[length(edges)+1]] <- c(j,i)
      } else if (rank_j - rank_i == 1){
        if (rank_i == 0) { edges[[length(edges)+1]] <- c(i,j)
        } else if (is_sub_flat(flat_i, flat_j)) edges[[length(edges)+1]] <- c(i,j)
      }
    }
  }
  # Build graph
  g <- igraph::graph_from_edgelist(el = do.call(rbind, edges), directed = TRUE)
  if (vcount(g) <= 10){
    # Optionally set nice labels
    V(g)$label <- sapply(flats, function(f) paste0("{", paste(f, collapse=","), "}"))
    # Plot using Sugiyama layout (hierarchical)
    plot(g, layout = layout_with_sugiyama(g)$layout, 
      vertex.size=8, vertex.label.cex=0.7,
      vertex.label.dist=2, vertex.label.degree=pi/2, edge.arrow.size = 0.3)
  } else{
    plot(g, layout = layout_with_sugiyama(g)$layout, vertex.size=8, edge.arrow.size=0.3, vertex.label=NA)
  }
  return(g)
}
G_MaxCones <- function(W, labels=NULL){
  # Maximum Cone Exchange
  n <- length(W)
  A <- matrix(0, n, n)
  for (i in 1:(n-1)){
    Wi <- W[[i]]
    for (j in (i+1):n){
      Wj <- W[[j]]
      d_ij <- sum(apply(Wi, 1, function(x) min(hamming(x, Wj))) > 0)==1
      A[i,j] <- d_ij
    }
  }
  g <- igraph::graph_from_adjacency_matrix(A + t(A), mode = "undirected")
  if (!is.null(labels) & vcount(g) <= 10) {
    V(g)$label <- sapply(labels, function(x) paste0("{", x, "}"))
    plot(g, vertex.size=3, vertex.label.cex=0.7, vertex.label.dist=1, vertex.label.degree=pi/2,)
  } else  plot(g, vertex.size=3, vertex.label=NA)
  return(g)
}

#### Fermat-Weber and Intersections ####

# Computes the FW objective function value of a point x for a sample V.
tr.sumDist <- function(V, x, symmetric=TRUE, to.x=TRUE){
  # Dependencies: Rfast
  # Note: If asymmetric, modifying "to.x" changes between max/min.
  if (class(V)[1] != "matrix") V <- Rfast::data.frame.to_matrix(V)
  if (class(x)[1] != "mumeric") x <- as.numeric(x)
  if (to.x) { vx <- Rfast::eachrow(-V, x, oper = "+")
  } else vx <- Rfast::eachrow(V, x, oper = "-")
  if (symmetric){
    vx.minmax <- Rfast::rowMinsMaxs(vx)
    return(sum(vx.minmax[2,]-vx.minmax[1,]))
  } else{
    n <- ncol(vx)
    return(sum(Rfast::rowsums(vx) - n*Rfast::rowMins(vx,value=TRUE)))
  }
}
# Solves the symmetric FW problem via min-cost flow LP.
tr.fw.mcf <- function(V){
  n <- nrow(V)
  d <- ncol(V)
  EL <- data.table(i = c(as.vector(row(V)) + d, as.vector(col(V))), 
    j = c(as.vector(col(V)), as.vector(row(V)) + n+d), 
    weight = c(as.vector(V), as.vector(-V)))
  G <- graph_from_data_frame(EL)
  edges <- as_edgelist(G, names = FALSE)
  IM <- matrix(0, nrow=vcount(G), ncol=ecount(G))
  for (j in seq_len(nrow(EL))){
    IM[edges[j,1],j] <- 1
    IM[edges[j,2],j] <- -1
  }
  rownames(IM) <- V(G)$name
  colnames(IM) <- paste0("e", seq_len(nrow(EL)))
  supply <- (rownames(IM) %in% seq(from=d+1,to=d+n)) - 
    (rownames(IM) %in% seq(from=d+n+1,to=d+2*n))
  mcf_LP <- lp("min", objective.in = E(G)$weight, const.mat = IM,
    const.dir = rep("==", nrow(IM)), const.rhs = supply, all.bin = TRUE)
  d0 <- -mcf_LP$objval
  w_flow <- which(mcf_LP$solution %==% 1)
  G_x <- igraph::reverse_edges(G, E(G)[w_flow])
  E(G_x)$weight[w_flow] <- -E(G)$weight[w_flow]
  x_1 <- which(V(G_x)$name == as.character(1))
  x_d <- which(V(G_x)$name %in% as.character(seq(2,d)))
  x_N <- which(V(G_x)$name %in% as.character(seq(d)))
  K_star <- tryCatch(igraph::distances(G_x, v = x_N, to = x_N, mode = "out"), error = function(e) NULL)
  if (is.null(K_star)) { warning("Degenerate LP solution.")
  } else {
    K_star <- K_star[order(as.integer(colnames(K_star))), order(as.integer(colnames(K_star)))]
    if (!all(apply(K_star, 1, function(x) tr.sumDist(V, x)) %<=% d0)) stop("BAD KSTAR")
    n_pi <- igraph::distances(G_x, v = x_1, to = x_d, mode = "out")
  }
  x0 <- c(0, n_pi[order(as.integer(colnames(n_pi)))])
  return(list(x0 = x0, d0 = d0, K_star = K_star))
}
# An alternative method for determining if a point is in FW(V) using gradients.
tr.fw.grad <- function(V, x){
  Vx <- Rfast::transpose(eachrow(V, x, oper = "-"))
  Vx.minmax <- Rfast::colMinsMaxs(Vx)
  T_max <- Rfast::eachrow(Vx, Vx.minmax[1,], oper = "-") %==% 0
  T_min <- Rfast::eachrow(Vx, Vx.minmax[2,], oper = "-") %==% 0
  T_colmax <- Rfast::colsums(T_max)
  T_colmin <- Rfast::colsums(T_min)
  if (max(c(T_colmax, T_colmin)) == 1){
    t_max <- Rfast::rowsums(T_max)
    t_min <- Rfast::rowsums(T_min)
    return(t_max - t_min)
  } else{
    n <- nrow(V); d <- ncol(V)
    SY <- data.table(i = 0, j = seq(d+1,d+n))
    ZT <- data.table(i = seq(d+n+1,d+2*n), j = d+2*n+1)
    T_max[!T_max] <- NA
    T_min[!T_min] <- NA
    YX <- data.table(reshape2::melt(T_max, na.rm = TRUE))
    XZ <- data.table(reshape2::melt(T_min, na.rm = TRUE))
    YX[, `:=` (i = Var2 + d, j = Var1)]
    XZ[, `:=` (i = Var1, j = Var2 + d+n)]
    EL <- rbindlist(list(SY, YX[,4:5], XZ[,4:5], ZT))
    EL[, capacity := 1]
    G <- igraph::graph_from_data_frame(EL)
    mf <- max_flow(G, source=V(G)["0"], target=V(G)[as.character(d+2*n+1)])
    if (mf$value == nrow(V)){
      return(rep(0, length(x)))
    } else{
      T_max <- Rfast::eachrow(T_max+0, T_colmax, oper = "/")
      T_min <- Rfast::eachrow(T_min+0, T_colmin, oper = "/")
      t_max <- Rfast::rowsums(T_max, na.rm = TRUE)
      t_min <- Rfast::rowsums(T_min, na.rm = TRUE)
      return(t_max - t_min)
    }
  }
}
# Backbone of Algorithm 1. Finds optimal r distance from boundary of cone (if it exists).
max_w_min <- function(K, V, C_matrix){
  n <- nrow(V)
  d <- ncol(V)
  VD <- data.table(i = as.vector(row(V)), j = as.vector(col(V)), v_ij = as.vector(V))
  setorder(VD, i, j)
  DC <- data.table::CJ(a = c(-1,1), k = c(0,1), i = seq(n), j = seq(d))
  DC[, `:=` (i = i + k*n, j = j + n*2)]
  DC[, `:=` (r = rep(seq(n*d*2), 2), c = i)]
  DC[a == -1, c := j]
  DC[k == 1, a := -a]
  b_vec <- c(-VD$v_ij, VD$v_ij)
  c_vec <- c(rep(1,n), rep(-1,n), rep(0,d))
  c_dir <- c(rep(">=", n*d*2))
  LP <- lp(direction = "min", objective.in = c_vec, const.dir = c_dir,
    const.rhs = b_vec, compute.sens = TRUE,
    dense.const = Rfast::data.frame.to_matrix(DC[,c("r","c","a")]))
  d0 <- LP$objval
  x0 <- tail(LP$solution, d)
  rho <- head(LP$duals, length(b_vec))
  SEC <- data.table(jj = seq(choose(d,2)), RcppAlgos::comboGeneral(d, 2))
  DC2 <- data.table::CJ(a = c(-1,0,1), k = c(0,1), i = seq(n), jj = seq(choose(d,2)))
  DC2[SEC, on = list(jj), `:=` (j1 = i.V1, j2 = i.V2)]
  DC2[k == 1, `:=` (j1 = j2, j2 = j1)]
  DC2[VD, on = list(i,j1=j), v_ij := i.v_ij]
  DC2[VD, on = list(i,j2=j), v_ik := i.v_ij]
  DC2[, `:=` (b = v_ij - v_ik, r = rep(seq(n*(d-1)*d), 3), c = i)]
  DC2[a == 1, c := j1 + n]
  DC2[a == -1, c := j2 + n]
  
  DC3 <- dcast(DC, r + i + j ~ c, value.var = "a", fill=0)[rho %>>% 0]
  DC3[, `:=` (k = as.numeric(i > n), jj = j-2*n)]
  DC3[k == 1, i := i-n]
  
  DC2[DC3[k==0], on = list(i,j1=jj), keep1 := TRUE]
  DC2[DC3[k==1], on = list(i,j2=jj), keep2 := TRUE]
  DC4 <- subset(DC2, keep1 == TRUE | keep2 == TRUE)
  DC4[a == 0, a := 1]
  
  A2_mat <- Matrix::sparseMatrix(i = DC4$r, j = DC4$c, x = DC4$a)[unique(DC4$r),]
  b2_vec <- DC4[a == -1]$b
  FW.hrep <- makeH(a1=as.matrix(A2_mat[,(n+1):(n+d)]), b1=-b2_vec)
  FW.hrep <- redundant(FW.hrep)$output
  
  K.hrep <- scdd(makeV(rays = -K))$output
  hrep <- rbind(K.hrep, FW.hrep)
  
  A_mat <- hrep[,-(1:2)]
  r <- as.numeric(rowMaxs(abs(A_mat), value=TRUE))
  r[which(K.hrep[,1] == 1)] <- 0
  r[which(seq(nrow(A_mat)) > nrow(K.hrep))] <- 0
  rA <-  cbind(r, A_mat)
  c_dir <- rep("<=", nrow(A_mat))
  c_dir[which(K.hrep[,1] == 1)] <- "=="
  c_dir[which(FW.hrep[,1] == 1) + nrow(K.hrep)] <- "=="
  b_vec <- c(K.hrep[,2], FW.hrep[,2])
  c_vec <- c(1, rep(0, ncol(K)))
  
  lb <- list(lower = list(ind = seq(ncol(rA)), val = c(0,rep(-Inf,ncol(K)))))
  LP <- Rglpk_solve_LP(obj = c_vec, mat = rA, dir = c_dir, rhs = b_vec, bounds = lb, max=TRUE)
  if (LP$status != 0) return(NULL)
  if (LP$optimum %<=% 0) return(NULL)
  x <- LP$solution[-1]
  if (!in_cone(x, K)) {
    SU <<- V
    x <<- x
    K <<- K
    stop("x not in K.")
  }
  if (!is_M_ultrametric(x, C_matrix)) stop("x not in B(M)")
  if (any(tr.fw.grad(V, x) %!=% 0)) stop("x not in FW(S)")
  if (w_min_C(x, C_matrix) %==% 0) return(NULL)
  return(x)
}
# Backbone of Algorithm 2. Finds point in intersection of pre-image(K) and FW(V).
does_fw_int_K <- function(V, V_matrix, C_matrix, ki, Ki, K_code, FWU_inK, preImage_Set){
  n <- nrow(V)
  d <- ncol(V)
  EL <- data.table(i = c(as.vector(row(V)) + d, as.vector(col(V))), 
    j = c(as.vector(col(V)), as.vector(row(V)) + n+d), 
    weight = c(as.vector(V), as.vector(-V)))
  G <- graph_from_data_frame(EL)
  IM <- get_incidence_matrix(G, directed = TRUE)
  supply <- (rownames(IM) %in% seq(from=d+1,to=d+n)) - 
    (rownames(IM) %in% seq(from=d+n+1,to=d+2*n))
  mcf_LP <- lp("min", objective.in = E(G)$weight, const.mat = IM,
    const.dir = rep("==", nrow(IM)), const.rhs = supply, all.bin = TRUE)
  d0_mcf <- -mcf_LP$objval
  w_flow <- which(mcf_LP$solution %==% 1)
  G_x <- igraph::reverse_edges(G, E(G)[w_flow])
  E(G_x)$weight[w_flow] <- -E(G)$weight[w_flow]
  x_1 <- which(V(G_x)$name == as.character(1))
  x_d <- which(V(G_x)$name %in% as.character(seq(2,d)))
  x_N <- which(V(G_x)$name %in% as.character(seq(d)))
  K_star <- tryCatch(igraph::distances(G_x, v = x_N, to = x_N, mode = "out"), error = function(e) NULL)
  if (is.null(K_star)) warning("Degenerate MCF LP solution.")
  if (!is.null(K_star)){
    K_star <- K_star[order(as.integer(colnames(K_star))), order(as.integer(colnames(K_star)))]
    n_pi <- igraph::distances(G_x, v = x_1, to = x_d, mode = "out")
    x0_mcf <- c(0, n_pi[order(as.integer(colnames(n_pi)))])
    trop_vert <- rbind(K_star, -t(K_star))
    if (all(apply(trop_vert, 1, function(x) tr.fw.grad(V, x)) %==% 0)){
      trop_vert_proj <- t(apply(trop_vert, 1, function(x) proj.tconv(x, V_matrix)))
      trop_vert_proj_codes <- apply(trop_vert_proj, 1, function(x) topological_profile(x, C_matrix))
      if (any(trop_vert_proj_codes == K_code)) return(TRUE)
      FWU_inK_proj <- proj.tconv(FWU_inK, -t(K_star))
      if (all(tr.fw.grad(V, FWU_inK_proj) %==% 0)){
        FWU_inK_projU <- proj.tconv(FWU_inK_proj, V_matrix)
        FWU_inK_projU_K <- topological_profile(FWU_inK_projU, C_matrix)
        if (FWU_inK_projU_K == K_code) return(TRUE)
      } else warning("Bad projection onto FW(S+Z).")
    } else warning("Bad K-Star (not FW points).")
  }
  VD <- data.table(i = as.vector(row(V)), j = as.vector(col(V)), v_ij = as.vector(V))
  setorder(VD, i, j)
  DC <- data.table::CJ(a = c(-1,1), k = c(0,1), i = seq(n), j = seq(d))
  DC[, `:=` (i = i + k*n, j = j + n*2)]
  DC[, `:=` (r = rep(seq(n*d*2), 2), c = i)]
  DC[a == -1, c := j]
  DC[k == 1, a := -a]
  A_mat <- Matrix::sparseMatrix(i = DC$r, j = DC$c, x = DC$a)
  b_vec <- c(-VD$v_ij, VD$v_ij)
  c_vec <- c(rep(1,n), rep(-1,n), rep(0,d))
  c_dir <- c(rep(">=", n*d*2))
  LP <- lp(direction = "min", objective.in = c_vec, const.dir = c_dir,
    const.rhs = b_vec, compute.sens = TRUE,
    dense.const = Rfast::data.frame.to_matrix(DC[,c("r","c","a")]))
  d0 <- LP$objval
  x0 <- tail(LP$solution, d)
  rho <- head(LP$duals, length(b_vec))
  
  SEC <- data.table(jj = seq(choose(d,2)), RcppAlgos::comboGeneral(d, 2))
  DC2 <- data.table::CJ(a = c(-1,0,1), k = c(0,1), i = seq(n), jj = seq(choose(d,2)))
  DC2[SEC, on = list(jj), `:=` (j1 = i.V1, j2 = i.V2)]
  DC2[k == 1, `:=` (j1 = j2, j2 = j1)]
  DC2[VD, on = list(i,j1=j), v_ij := i.v_ij]
  DC2[VD, on = list(i,j2=j), v_ik := i.v_ij]
  DC2[, `:=` (b = v_ij - v_ik, r = rep(seq(n*(d-1)*d), 3), c = i)]
  DC2[a == 1, c := j1 + n]
  DC2[a == -1, c := j2 + n]
  DC3 <- dcast(DC, r + i + j ~ c, value.var = "a", fill=0)[rho %>>% 0]
  DC3[, `:=` (k = as.numeric(i > n), jj = j-2*n)]
  DC3[k == 1, i := i-n]
  DC2[DC3[k==0], on = list(i,j1=jj), keep1 := TRUE]
  DC2[DC3[k==1], on = list(i,j2=jj), keep2 := TRUE]
  DC4 <- subset(DC2, keep1 == TRUE | keep2 == TRUE)
  DC4[a == 0, a := 1]
  A2_mat <- Matrix::sparseMatrix(i = DC4$r, j = DC4$c, x = DC4$a)[unique(DC4$r),]
  b2_vec <- DC4[a == -1]$b
  FW.hrep <- makeH(a1=as.matrix(A2_mat[,(n+1):(n+d)]), b1=-b2_vec)
  FW.rdnd <- redundant(FW.hrep)$output
  FW.scdd <- tryCatch(scdd(FW.rdnd$output)$output, error = function(e) NULL)
  if (!is.null(FW.scdd)){
    FW.set <- FW.scdd[which(FW.scdd[,1] == 0),-(1:2), drop = FALSE]
    if (nrow(FW.set) > 0) {
      FW.set_proj <- t(apply(FW.set, 1, function(x) proj.tconv(x, V_matrix)))
      FW.set_proj_code <- apply(FW.set_proj, 1, function(x) topological_profile(x, C_matrix))
      if (any(FW.set_proj_code == K_code)) return(TRUE)
    }
  } # else warning("Error with SCDD.") # Silencing for now.
  for (Hspace in preImage_Set){
    A_mat <- rbind(FW.rdnd[,-(1:2)], -Hspace)
    r <- as.numeric(rowMaxs(abs(A_mat), value=TRUE))
    r[which(FW.rdnd[,1] == 1)] <- 0
    #r[which(seq(nrow(A_mat)) > nrow(FW))] <- 0
    rA <-  cbind(r, A_mat)
    c_dir <- rep("<=", nrow(A_mat))
    c_dir[which(FW.rdnd[,1] == 1)] <- "=="
    b_vec <- c(FW.rdnd[,2], rep(0, nrow(Hspace)))
    c_vec <- c(1, rep(0, ncol(V)))
    lb <- list(lower = list(ind = seq(ncol(rA)), val = c(0,rep(-Inf,ncol(SZ)))))
    LP <- Rglpk_solve_LP(obj = c_vec, mat = rA, dir = c_dir, rhs = b_vec, bounds = lb, max=TRUE)
    if (LP$status != 0) next
    if (LP$optimum %==% 0) next
    x <- LP$solution[-1]
    x <- x-min(x)
    if (all(tr.fw.grad(V, x) %==% 0)){
      if (in_Hspace(x, Hspace)) { return(TRUE)
      } else warning("Bad intersection. Not in Preimage.")
    } else warning("Bad intersection. Not FW point.")
  }
  return(FALSE)
}
