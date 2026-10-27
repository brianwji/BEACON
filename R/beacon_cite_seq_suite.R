
#################################
#### Learn BEACON parameters ####
#################################

learn_beacon_params_cite <- function(data,thetas_b0,sig_max,min_count,max_iter,nb_params){

  #Initialize parameters
  m <- nrow(data)
  n <- ncol(data)
  counts_total <- colSums(data)

  #Initial parameters
  print("Solving initial parameters...")
  thetas_b <- thetas_b0
  output <- mle_alphas_thetas(data,thetas_b,sig_max,min_count,nb_params)
  thetas_s <- output$thetas_s_max
  LLs <- output$LLs_max
  LL_joint <- sum(LLs,na.rm = TRUE)
  alphas <- output$alphas_max
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
    output <- mle_alphas_thetas(data,thetas_b_test,sig_max,min_count,nb_params)
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
      cat(sprintf("Joint log-likelihood = %.0f\n", tail(LLs_joint, 1)))

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

mle_alphas_thetas <- function(data,thetas_b,sig_max,min_count,nb_params){

  #Size
  m <- nrow(data)
  n <- ncol(data)

  #Assess the probability of data coming from the background
  pdata <- matrix(0,nrow=m,ncol=n)

  for(i in 1:m){

    #Marginal
    x_i <- data[i,]

    #Probability of data given background
    p_i <- pnbinom(x_i-1, size = nb_params$size, mu = nb_params$mu,lower.tail = FALSE)
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
    pcut <- 0.2
    pinds <- which(pvec <= pcut)

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
