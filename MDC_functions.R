# discretization
discretize_corr = function(C,D){
  # Ensure C and D are square matrices of the same dimension
  if (!all(dim(C) == dim(D))) {
    stop("C and D must be square matrices of the same dimensions.")
  }
  
  # Initialize C_filter with zeros
  C_filter <- matrix(0, nrow = dim(C)[1], ncol = dim(C)[2])
  
  # Identify indices of the upper triangular off-diagonal part
  upper_tri_indices <- which(upper.tri(C), arr.ind = TRUE)
  
  # Vectorized operations for efficiency
  positive_indices <- which(C[upper_tri_indices] > 0 & D[upper_tri_indices] == 1)
  negative_indices <- which(C[upper_tri_indices] < 0 & D[upper_tri_indices] == 1)
  na_indices <- which(is.na(C[upper_tri_indices]) )
  
  # Update C_filter based on conditions
  # print(upper_tri_indices[positive_indices])
  C_filter[upper_tri_indices[positive_indices,]] <- 1
  C_filter[upper_tri_indices[negative_indices,]] <- -1
  C_filter[upper_tri_indices[na_indices,]] <- NA
  
  return(C_filter)
}

# Calculate P1 transition matrix
transition_P1 = function(corr_mat1,corr_mat2){
  P_m <- matrix(0, nrow = 3, ncol = 3)
  for (i in 1:(nrow(corr_mat2)-1)) {
    for (j in (i+1):ncol(corr_mat2)) {
      M1_val <- corr_mat1[i, j]
      M2_val <- corr_mat2[i, j]
      # Increment the corresponding element in P_m
      # print(i)
      # print(j)
      if (!is.na(M2_val)) {
        if (M2_val == -1) {
          P_m[M1_val + 2, 1] <- P_m[M1_val + 2, 1] + 1
          } 
        else if (M2_val == 0) {
          P_m[M1_val + 2, 2] <- P_m[M1_val + 2, 2] + 1
          } 
        else if (M2_val == 1) {
          P_m[M1_val + 2, 3] <- P_m[M1_val + 2, 3] + 1
          }
        }
    }
  }
  # Normalize each row so that the sum of each row equals 1
  # P_m <- sweep(P_m, 1, rowSums(P_m), "/")
  # P_m[is.na(P_m)] <- 0
  return(P_m)
}

# Calculate P1 transition matrix
transition_P2 = function(corr_mat1,corr_mat2,corr_mat3){
  P_m <- matrix(0, nrow = 3^2, ncol = 3)
  for (i in 1:(nrow(corr_mat2)-1)) {
    for (j in (i+1):ncol(corr_mat2)) {
      M1_val <- corr_mat1[i, j]
      M2_val <- corr_mat2[i, j]
      M3_val <- corr_mat3[i, j]
      # Increment the corresponding element in P_m
      # ifelse(M3_val==-1, m3=1, ifelse(M3_val==0, m3=2, m3=3))
      if (!is.na(M1_val) & !is.na(M2_val)){
        if (M1_val==-1) {m1=1} else if (M1_val==0) {m1=2} else {m1=3}
        if (M2_val == -1) {
          P_m[m1*3-2, M3_val + 2] <- P_m[m1*3-2, M3_val + 2] + 1
          } 
        else if (M2_val == 0) {
          P_m[m1*3-1, M3_val + 2] <- P_m[m1*3-1, M3_val + 2] + 1
          } 
        else if (M2_val == 1) {
          P_m[m1*3,   M3_val + 2] <- P_m[m1*3,   M3_val + 2] + 1
          }
        }
    }
  }
  # # Normalize each row so that the sum of each row equals 1
  # P_m <- sweep(P_m, 1, rowSums(P_m), "/")
  # P_m[is.na(P_m)] <- 0
  return(P_m)
}

# Calculate Fisher Z-score
Fisher_Z = function(r){
  Z=log((1+r)/(1-r))/2
  return(Z)
}

# Test correlation difference significance by Fisher Z-score
Test_FisherZ = function(Z1,Z2,n1,n2){
  z=(Z1-Z2)/(1/(n1-3)+1/(n2-3))
  return(z)
}
