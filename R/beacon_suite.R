#' Demultiplexing Single-cell RNA-sequencing Experiments Using BEACON
#'
#' Performs sample assignment by estimating and removing background barcode oligos
#'
#'
#' @param data A numeric matrix containing HTO counts,
#'   with barcodes as rows and cells as columns.
#' @param params A list of BEACON parameters.
#'
#' @return A named numeric vector of sample assignments.
#'
#' @export

######################
#### BEACON calls ####
######################

beacon_calls <- function(data,params){

  #Parameters
  sig_max <- params$sig_max
  min_count <- params$min_count
  max_bf <- params$max_bf
  max_iter <- params$max_iter

  #Initialize background distribution
  m <- nrow(data)
  n <- ncol(data)
  thetas_b0 <- rowSums(data) / sum(rowSums(data))

  #Estimate the doublet fraction
  dprior <- doublet_prior(n)

  #Learn beacon parameters
  output <- learn_beacon_params(data,thetas_b0,sig_max,min_count,max_iter)

  #Remove the background
  thetas_b <- output$thetas_b
  alphas <- output$alphas
  thetas_b_mat <- matrix(thetas_b,m,n,byrow=FALSE)
  alphas_mat <- matrix(alphas,m,n,byrow=TRUE)
  counts <- colSums(data)
  counts_mat <- matrix(counts,m,n,byrow=TRUE)
  thetas_s <- output$thetas_s
  data_sig <- ceiling(counts_mat * (1-alphas_mat) * thetas_s)

  #Assign labels

  #Second most abundant barcode
  top2 <- apply(data_sig, 2, function(x)
    sort(x, decreasing = TRUE)[2]
  )

  #Sort the second most abundant barcode
  top2_order <- order(top2,decreasing=TRUE)

  #Call doublets
  n_doublets <- ceiling (dprior * n)
  doublet_inds <- top2_order[1:n_doublets]

  #Negatives
  alpha_cut <- max_bf
  sigs <- output$sigs
  neg_inds <- which(sigs == 0)
  neg_inds <- unique(c(neg_inds,which(alphas > alpha_cut)))

  #Call singlets
  top1_inds <- apply(data_sig, 2, function(x)
    order(x, decreasing = TRUE)[1]
  )

  #Final calls
  labels_beacon <- top1_inds
  labels_beacon[doublet_inds] <- 1000
  labels_beacon[neg_inds] <- 0
  names(labels_beacon) <- colnames(data)

  return(labels_beacon)

}





#################################
#### Learn BEACON parameters ####
#################################

learn_beacon_params <- function(data,thetas_b0,sig_max,min_count,max_iter){

  #Initialize parameters
  m <- nrow(data)
  n <- ncol(data)
  counts_total <- colSums(data)

  #Initial parameters
  print("Solving initial parameters...")
  thetas_b <- thetas_b0
  output <- mle_alphas_thetas(data,thetas_b,sig_max,min_count)
  thetas_s <- output$thetas_s_max
  LLs <- output$LLs_max
  LL_joint <- sum(LLs,na.rm = TRUE)
  alphas <- output$alphas_max
  sigs <- output$sigs
  del <- 100
  count <- 0
  tol = 1 #For convergence

  #Storing
  LLs_joint <- LL_joint
  thetas_b_iter <- thetas_b #keep track of iterations
  thetas_s_iter <- list() #keep track of iterations
  thetas_s_iter[[1]] <- thetas_s
  alphas_iter <- t(alphas)

  while(del > tol){

    count <- count + 1

    #Optimize background count distribution
    alpha_mat <- matrix(alphas,m,n,byrow=TRUE)
    count_mat <- matrix(counts_total,m,n,byrow=TRUE)
    thetas_b_test_un <- rowMeans(data - count_mat * (1 - alpha_mat) * thetas_s, na.rm=TRUE)

      #No negative values
      thetas_b_test_un[thetas_b_test_un < 0] = 0
      thetas_b_test <- thetas_b_test_un / sum(thetas_b_test_un)

    #Optimize thetas_s and alphas
    output <- mle_alphas_thetas(data,thetas_b_test,sig_max,min_count)
    thetas_s_test <- output$thetas_s_max
    LLs_test <- output$LLs_max
    LL_joint_test <- sum(LLs_test,na.rm = TRUE)
    LLs_joint_test <- c(LLs_joint,LL_joint_test)
    alphas_test <- output$alphas_max
    sigs_test <- output$sigs

    #Calculate difference in the likelihood function
    diffLL_test <- diff(LLs_joint_test)
    del <- diffLL_test[length(diffLL_test)]

    #Accept criteria
    if(del > tol){

      print(paste0("Iteration ", count, "..."))

      #Update signal
      thetas_s <- thetas_s_test
      thetas_s_iter[[count+1]] <- thetas_s

      #Update background
      thetas_b <- thetas_b_test
      thetas_b_iter <- cbind(thetas_b_iter,thetas_b)

      #Update alphas
      alphas <- alphas_test
      alphas_iter <- cbind(alphas_iter,t(alphas))

      #Update log likelihood
      LLs_joint <- LLs_joint_test
      cat(sprintf("Joint log-likelihood = %.0f\n", utils::tail(LLs_joint, 1)))

      #Update sigs
      sigs <- sigs_test

    }

    #Maximum iteration
    if(count > max_iter){
      print(paste("Maximum iterations of",max_iter,"reached"))
      break
    }

  }

  #Convergence status
  if(count < max_iter & del < tol ){
    print(paste("Parameters converged!"))
  }

  output <- list()
  output$thetas_s <- thetas_s
  output$alphas <- alphas
  output$thetas_b <- thetas_b
  output$thetas_b_iter <- thetas_b_iter
  output$LLs_joint <- LLs_joint
  output$sigs <- sigs
  output$thetas_s_iter <- thetas_s_iter
  output$alphas_iter <- alphas_iter


  return(output)
}



#################################
#### MLE alphas and thetas_s ####
#################################

mle_alphas_thetas <- function(data,thetas_b,sig_max,min_count){

  #Size
  m <- nrow(data)
  n <- ncol(data)

  #Normalize log abundances
  logdata <- log(data+1)
  logdata <- logdata - matrix(colMeans(logdata),nrow=m,ncol=n,byrow=TRUE)

  #EM to learn the background distribution
  means_em <- matrix(0,nrow=m,ncol=2)
  stds_em <- matrix(0,nrow=m,ncol=2)
  mix_em <- matrix(0,nrow=m,ncol=2)
  pdata <- matrix(0,nrow=m,ncol=n)

  for(i in 1:m){

    #Marginal
    x_i <- logdata[i,]

    #Initialize mixture proportions
    mix_low <- 1 - 1/m
    mix_high <- 1/m

    #Initialize means and variances
    sort_inds <- order(x_i)
    inds_low <- sort_inds[1:floor(mix_low * n)]
    inds_high <- sort_inds[ceiling(mix_low * n):n]
    low_i <- x_i[inds_low]
    high_i <- x_i[inds_high]
    mean_low <- mean(low_i)
    mean_high <- mean(high_i)
    std_low <- stats::sd(low_i)
    std_high <- stats::sd(high_i)

    err = 1
    tol = 1e-5 #For convergence
    while(err > tol){

      #Expectation
      w_low = mix_low * stats::dnorm(x_i, mean = mean_low, sd = std_low)
      w_high = mix_high * stats::dnorm(x_i, mean = mean_high, sd = std_high)
      w_low_norm <- w_low / (w_low + w_high)
      w_high_norm <- w_high / (w_low + w_high)

      #Dirichlet prior on mixture fractions
      alpha_total <- 500
      alpha_high <- alpha_total / m
      alpha_low <- alpha_total - alpha_high

      #Maximization
      N_low <- sum(w_low_norm)
      N_high <- sum(w_high_norm)
      mix_low_new <- (N_low + alpha_low  -1) / (N_low + N_high + alpha_low + alpha_high -2)
      mix_high_new <- (N_high + alpha_high - 1) / (N_low + N_high + alpha_low + alpha_high - 2)
      mean_low_new <- 1 / (n * mix_low_new) * sum(w_low_norm * x_i)
      mean_high_new <- 1 / (n * mix_high_new) * sum(w_high_norm * x_i)
      var_low_new <- 1 / (n * mix_low_new) * sum(w_low_norm * (x_i - mean_low_new)^2)
      var_high_new <- 1 / (n * mix_high_new) * sum(w_high_norm * (x_i - mean_high_new)^2)
      std_low_new <- sqrt(var_low_new)
      std_high_new <- sqrt(var_high_new)

      #Calculate error
      e1 <- c(mean_low, mean_high, std_low, std_high)
      e2 <- c(mean_low_new, mean_high_new, std_low_new, std_high_new)
      err = sqrt(mean((e2-e1)^2))

      #Update
      mix_low <- mix_low_new
      mix_high <- mix_high_new
      mean_low <- mean_low_new
      mean_high <- mean_high_new
      std_low <- std_low_new
      std_high <- std_high_new

    }

    #Storing
    mix_em[i,] <- c(mix_low, mix_high)
    means_em[i,] <- c(mean_low, mean_high)
    stds_em[i,] <- c(std_low, std_high)

    p_i <- stats::pnorm(x_i, mean = mean_low, sd = std_low, lower.tail = FALSE)
    pdata[i,] <- p_i


  }

  #Find the optimal set of alphas and thetas_s
  alpha_scan <- seq(0,0.998,by = 0.002)
  thetas_s_max <- matrix(0, nrow = m, ncol = n)
  alphas_max <- matrix(0, nrow = 1,ncol = n)
  LLs_max <- matrix(0, nrow = 1, ncol = n)
  sigs <- matrix(0,nrow = 1, ncol = n)

  for(i in 1:n){

    #Identify the candidate signal barcodes
    vec <- data[,i]
    pvec <- pdata[,i]

    #Probability cut off
    pcut <- 0.05
    pinds <- which(pvec < pcut)

    #Count cut off
    count_inds <- which(vec >= min_count)

    #Take the intersection
    ainds <- intersect(pinds,count_inds)

    asc_inds <- order(pvec[ainds])
    k <- length(ainds)
    sigs[i] <- k

    #If candidate signal barcodes < maximum allowable
    if(k < sig_max){
      sig_m <- k
      inds_s <- ainds[asc_inds[1:sig_m]]
      inds_b <- setdiff(1:m,inds_s)
    }

    #If candidate signal barcodes >= maximum allowable
    if(k >= sig_max){
      sig_m <- sig_max
      inds_s <- ainds[asc_inds[1:sig_m]]
      inds_b <- setdiff(1:m,inds_s)
    }

    #If no candidate signal barcodes
    if(k == 0){
      alphas_max[i] <- 1
      next
    }

    #Inputs into the multinomial mixture model
    params <- list()
    params$alpha_scan <- alpha_scan
    params$thetas_b <- thetas_b
    params$vec <- vec
    params$inds_s <- inds_s
    params$sig_m <-sig_m

    #Find the optimal parameters by scanning over alphas
    X <- mle_alpha_theta(params)

    #Store
    thetas_s_max[inds_s,i] <- X$theta_s_max
    alphas_max[i] <- X$alpha_max
    LLs_max[i] <- X$LL_max

  }

  #No negative thetas_s
  thetas_s_max[thetas_s_max < 0] <- 0
  thetas_s_max <- thetas_s_max / matrix(colSums(thetas_s_max),m,ncol(thetas_s_max),byrow=TRUE)

  #Output
  output <- list()
  output$thetas_s_max <- thetas_s_max
  output$alphas_max <- alphas_max
  output$LLs_max <- LLs_max
  output$sigs <- sigs

  return(output)

}





###############################
#### MLE alpha and theta_s ####
###############################

mle_alpha_theta <- function(params){

  #Initial parameters
  alpha_scan <- params$alpha_scan
  n <- length(alpha_scan)
  thetas_b <- params$thetas_b
  vec <- params$vec
  inds_s <- params$inds_s
  sig_m <- params$sig_m
  x <- vec[inds_s]
  a <- thetas_b[inds_s]
  A <- sum(a)
  X <- sum(x)
  m <- length(vec)

  #Optimal thetas_s as a function of alpha_scan
  xs <- matrix(x,sig_m,n,byrow=FALSE)
  as <- matrix(a,sig_m,n,byrow=FALSE)
  alphas <- matrix(alpha_scan,sig_m,n,byrow=TRUE)
  thetas_scan <- (xs * (1 - alphas + alphas * A) - alphas * X * as ) / (X * (1-alphas))

  #Calculate the likelihood of each set of alpha and theta_s
  alphas_mat <- matrix(alpha_scan,m,n,byrow=TRUE)
  thetas_b_mat <- matrix(thetas_b,m,n,byrow=FALSE)
  thetas_s_mat <- matrix(0,m,n)
  thetas_s_mat[inds_s,] <- thetas_scan
  p <- alphas_mat * thetas_b_mat + (1-alphas_mat) * thetas_s_mat
  XX <- matrix(vec,m,n,byrow=FALSE)
  LLs <- mnloglikelihood(XX,p)

  #Find the set of alphas and thetas_s that maximize the log likelihood
  max_ind <- which.max(LLs)
  if(length(max_ind)==0){
    max_ind = 1
  }
  alpha_max <- alpha_scan[max_ind]
  theta_s_max <- thetas_scan[,max_ind]
  LL_max <- LLs[max_ind]


  #Output
  output <- list()
  output$alpha_max <- alpha_max
  output$theta_s_max <- theta_s_max
  output$LL_max <- LL_max
  return(output)


}

#####################################
#### Multinomial log-likelihood #####
#####################################

mnloglikelihood <- function(XX,p){
 m <- nrow(XX)
 n <- ncol(XX)
 count <- sum(XX[,1])

 #Term 1
 term1 <- sumlog(count)
 terms1  <- rep(term1,times=n)

 #Term 2
 terms2 <- colSums(XX * log(p))

 #Term 3
 x1 <- XX[,1]
 term3 <- 0
 for(j in 1:m){
 term3 <- term3 + sumlog(max(1,x1[j]))
 }
 terms3 <- rep(term3,times=n)

 #Log-likelihood
 LLs <- terms1 + terms2 - terms3
 return(LLs)
}


#################
#### Sum log ####
#################

sumlog <- function(x){
  Y <- sum(log(1:x))
  return(Y)
}

###################################
#### Doublet prior estimation #####
###################################

doublet_prior <- function(n_cells){

  #Taken from 10X v3 user guide
  dps <- c(0.004,0.008,0.016,0.023,0.031,0.039,0.046,0.054,0.061,0.069,0.076)
  cells_recovered <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)

  #Regress
  x <- as.vector(cells_recovered)
  y <- as.vector(dps)
  fit <- stats::lm(y ~ x)

  # Predict
  dprior <- stats::predict(fit, newdata = data.frame(x = n_cells))
  return(dprior)
}
