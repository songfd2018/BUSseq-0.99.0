############################################################
# Internal functions to be used in the BUSseq_MCMC function#
############################################################
#get the posterior mode
.getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Calculate the BIC function of BUSseq
.cal_BIC_BUSseq<-function(BUSseqfits_obj){
  
  if(is(BUSseqfits_obj,"BUSseqfits")){
    Read_list <- BUSseqfits_obj$CountData_raw
    .B <- BUSseqfits_obj$n.batch
    .nb <- BUSseqfits_obj$n.perbatch
    .N <- BUSseqfits_obj$n.cell
    .G <- BUSseqfits_obj$n.gene
    .K <- BUSseqfits_obj$n.celltype
    .gamma <- BUSseqfits_obj$gamma.est
    .logmu <- BUSseqfits_obj$alpha.est + BUSseqfits_obj$beta.est
    .nu <- BUSseqfits_obj$nu.est
    .delta <- BUSseqfits_obj$delta.est
    .phi <- BUSseqfits_obj$phi.est
    .pi <- BUSseqfits_obj$pi.est
    
    Read <- NULL
    delta_vec <- NULL
    for(b in 1:.B){
      Read <- rbind(Read,t(Read_list[[b]]))
      delta_vec <- c(delta_vec,.delta[[b]])
    }
    
    args=list(Reads=as.integer(Read), K = as.integer(.K),  BNum = as.integer(.B), CNum = as.integer(.nb), GNum=as.integer(.G), #observed data
              mu = .logmu, nu= .nu, delta = delta_vec, gamma = .gamma, phi = .phi, pi = .pi) #parameter
    #largest.x = as.integer(largest_x)) 
    
    loglikelihood<-.Call('Cal_LogLike', args)
    bic <- (-2) * loglikelihood + log(.G * .N) * ((.B + .G) * .K + 2 * .B + .G * (.B * 2 - 1) + .N - .B)
    return(bic)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
  

}

#####################################
#Estimate intrinsic gene indicators #
#####################################
.fdrDEindicator <- function(PPI, kappa){
  .G <- nrow(PPI)
  .K <- ncol(PPI)
  
  .xi <- 1 - PPI[,2:.K]
  ind_intr <- .xi <= kappa
  fdr <- sum(.xi[ind_intr])/sum(ind_intr)
  
  return(fdr)
}

#calculate the DE posterior probability threshold
.postprob_DE_thr_fun <- function(PPI, fdr_threshold=0.1){
  
  .K <- ncol(PPI)
  kappa_fdr_matr <- NULL
  kappa_set <- rev(1 - sort(unique(c(PPI[,2:.K]))))
  for(kappa in kappa_set){
    fdr <- .fdrDEindicator(PPI, kappa=kappa)
    kappa_fdr_matr <- rbind(kappa_fdr_matr, c(kappa, fdr))
    if(fdr > fdr_threshold){
      break
    }
  }
  ind <- which(kappa_fdr_matr[ ,2] <= fdr_threshold)
  ind2 <- which.max(kappa_fdr_matr[ind,2])
  ind3 <- ind[ind2] # the index that has the maximum fdr
  # but less than fdr_threshold
  #message(c("Posterior probability threshold = ",
  #          1-as.numeric(kappa_fdr_matr[ind3,1]),"\n"))
  #message("The output is a scalar.\n")
  return(1-as.numeric(kappa_fdr_matr[ind3,1])) #1-kappa
}

#Estimate intrinsic gene indicators
.estimate_IG_indicators <- function(PPI, postprob_DE_threshold = 0.5){

  K<-ncol(PPI)
  EstL <- PPI[,2:K]
  EstL[PPI[,2:K] >= postprob_DE_threshold] <- 1
  EstL[PPI[,2:K] < postprob_DE_threshold] <- 0
  message("The output format is a matrix.\n")
  message(paste0("Each row represents a gene, and each column",
                 " corresponds to a cell type from 2 to K\n"))
  return(EstL)
}

#Intrinsic gene index
.IG_index <- function(EstIGindicators){
  ind <- which(rowSums(EstIGindicators) > 0)
  message(c(length(ind), " intrinsic genes are found.\n"))
  message("The output format is a vector implying the intrinsic gene",
          " indices.\n")
  return(ind)
}

##############################
# The MCMC Sampling Function #
##############################
BUSseq_MCMC <- function(ObservedData, n.celltypes, n.iterations = 2000, 
                     n.burnin = floor(n.iterations / 2),
                     n.unchanged = min(floor(n.iterations * 0.3), 500),
                     working_dir = ".", hyper_pi = NA,
                     hyper_mu_alpha = NA, hyper_mu_nu = NA,
                     hyper_sd_alpha = NA, hyper_sd_nu = NA,
                     hyper_sd_delta = NA, hyper_var_gamma = NA,
                     hyper_phi = NA, hyper_tau0 = NA, 
                     hyper_p = NA, hyper_slab = NA,
                     seed = 123,
                     showIteration = TRUE){
  
  Read <- NULL #input data Y
  nb <- NULL
  
if(is(ObservedData, "list")){         #The input data format is list
    #Each row represents a gene and each column is a cell
    B <- length(ObservedData)              #batch number    
    for(b in seq_len(B)){           #input data
      Read <- rbind(Read,t(ObservedData[[b]]))
      nb <- c(nb,ncol(ObservedData[[b]]))
    }
    ###################################
    ###Test the consistency of genes###
  }else{
    stop(paste0("ObservedData must be a \"list\" object!\n"))
  }
  
  if(B < 2){
    stop("The batch number must be greater than one.\n")
  }
  
  
  K <- n.celltypes
  G<-ncol(Read)
  N <- sum(nb)
  
  if(sum(K > nb) > 0){
    stop(paste0("The sample size in any batch must be greater",
                " than the assumed cell type number.\n"))
  }
  
  ####Record the posterior samplingo on the hard disk
  setwd(working_dir)
  sampling_dir <- paste0(working_dir,"/MCMC_sampling",sep="")
  dir.create(sampling_dir, showWarnings = F)
  
  #save.image(paste0(sampling_dir,"/raw_data.Rdata",sep=""))
  
  #create file for storing poserterior sampling
  alpha.file<-paste0(sampling_dir,"/alpha.txt",sep="")
  file.create(alpha.file)
  
  beta.file<-paste0(sampling_dir,"/beta.txt",sep="")
  file.create(beta.file)
  
  nu.file<-paste(sampling_dir,"/nu.txt",sep="")
  file.create(nu.file)
  
  delta.file<-paste(sampling_dir,"/delta.txt",sep="")
  file.create(delta.file)
  
  gamma.file<-paste(sampling_dir,"/gamma.txt",sep="")
  file.create(gamma.file)
  
  phi.file<-paste(sampling_dir,"/phi.txt",sep="")
  file.create(phi.file)
  
  pi.file<-paste(sampling_dir,"/pi.txt",sep="")
  file.create(pi.file)
  
  w.file <-paste(sampling_dir,"/w.txt",sep="")
  file.create(w.file)
  
  prob.file <- paste(sampling_dir,"/p.txt",sep="")
  file.create(prob.file)
  
  tau0.file <- paste(sampling_dir,"/tau0.txt",sep="")
  file.create(tau0.file)
  
  l.file <- paste(sampling_dir,"/l.txt",sep="")
  file.create(l.file)
  
  x.file <-paste(sampling_dir,"/x.txt",sep="")
  file.create(x.file)
  
  set.seed(seed)
  
  #########################
  #specify hyperparameters#
  #########################
  if(is.na(hyper_pi)){
    hyper_pi <- 2
  }
  if(is.na(hyper_sd_alpha)){
    hyper_sd_alpha <- sqrt(5)
  }
  if(is.na(hyper_sd_nu)){
    hyper_sd_nu <- sqrt(5)
  }
  if(is.na(hyper_sd_delta)){
    hyper_sd_delta <- sqrt(5)
  }
  if(any(is.na(hyper_var_gamma))){
    hyper_var_gamma <- c(3,3)
  }
  if(any(is.na(hyper_phi))){
    hyper_phi <- c(1,0.1)
  }
  if(any(is.na(hyper_p))){
    hyper_p <- c(1,3)
  }
  if(any(is.na(hyper_tau0))){
    hyper_tau0 <- c(1,0.1)
  }
  if(is.na(hyper_slab)){
    hyper_slab <- 100
  }
  
  xi <- hyper_pi
  sigma.a <- hyper_sd_alpha
  sigma.c <- hyper_sd_nu
  sigma.d <- hyper_sd_delta
  sigma.z <- hyper_var_gamma
  kappa <- hyper_phi[1]
  tau <- hyper_phi[2]
  a.tau <- hyper_tau0[1]
  b.tau <- hyper_tau0[2]
  a.p <- hyper_p[1]
  b.p <- hyper_p[2]
  tau1 <- hyper_slab#tau_1^2
  
  ####################
  #set initial values#
  ####################
  p.init <- runif(1,0,0.5)
  tau0.init <- 0.005
  gamma.init <- matrix(0,B,2)
  gamma.init[,2] <- 0.1
  pi.init <- matrix((K:1)/sum(1:K),K,B)
  phi.init <- matrix(5,B,G)
  
  ###Initial cell type indicators w_bi
  w.init <- rep(0,N)
  for(b in 1:B){
    for(i in 1:nb[b]){
      if(b==1){
        row.index <- i
      }else{
        row.index <- i + sum(nb[1:(b-1)])
      }
      w.init[row.index] <- which(rmultinom(1,1,pi.init[,b])==1) - 1
    }
  }
  
  ###Set data-dependent initial values
  delta.init <- matrix(0,N)
  sum_per_cell<-apply(Read,1,sum)
  for(b in 1:B){
    for(i in 2:nb[b]){
      if(b==1){
        row.index<-i
        delta.init[row.index] <- log(sum_per_cell[row.index]) - log(sum_per_cell[1])
      }else{
        row.index<-i + sum(nb[1:(b-1)])
        delta.init[row.index] <- log(sum_per_cell[row.index]) - log(sum_per_cell[1+sum(nb[1:(b-1)])])
      }
    }
  }
  
  raw_Means <- array(NA,dim=c(B,K,G))
  Reads_corrected<-matrix(NA,N,G)
  for(i in 1:N){
    Reads_corrected[i,]<-Read[i,]/exp(delta.init[i])
  }
  for(b in 1:B){
    if(b==1){
      index_batchb <- 1:nb[1]
    }else{
      index_batchb <- sum(nb[1:(b-1)]) + 1:nb[b]
    }
    read_b <- Reads_corrected[index_batchb,]
    for(k in 1:K){
      index_celltypek <- which(w.init[index_batchb]==k - 1)
      if(length(index_celltypek)>1){
        raw_Means[b,k,] <- apply(log(1 + read_b[index_celltypek,]),2,mean)
      }else if(length(index_celltypek)==1){
        raw_Means[b,k,] <- log(1 + read_b[index_celltypek,] )
      }else{
        raw_Means[b,k,] <- apply(log(1 + read_b),2,mean)
      }
    }
  }
  
  alpha.init <- raw_Means[1,1,]
  
  beta.init <- matrix(0,G,K)
  for(k in 2:K){
    beta.init[,k] <- raw_Means[1,k,] - alpha.init
  }
  
  l.init <- matrix(0,G,K)
  for(j in 1:G){
    for(k in 2:K){
      log_rat <- 0
      log_rat <- log_rat + log(p.init/(1-p.init))
      log_rat <- log_rat - log(tau1)/2 + log(tau0.init)/2
      log_rat <- log_rat - beta.init[j,k]^2/2/tau1 + beta.init[j,k]^2/2/tau0.init
      l.init[j,k] <- rbinom(1,1,1/(1+exp(-log_rat)))
    }
  }
  
  nu.init <- matrix(0,B,G)
  for(b in 2:B){
    nu.init[b,]<-apply(raw_Means[b,,],2,mean)-apply(raw_Means[1,,],2,mean)
  }
  
  
  z.init <- matrix(1,N,G)
  x.init <- Read
  
  ###Set the empirical estiamtors as the prior for alpha and nu
  if(any(is.na(hyper_mu_alpha))){
    mu.a <- alpha.init
  }else{
    mu.a <- hyper_mu_alpha
  }
  if(any(is.na(hyper_mu_nu))){
    mu.c<-apply(nu.init,1,mean)
  }else{
    mu.c <- hyper_mu_nu
  }
  
  ###############
  #MCMC sampling#
  ###############
  ###Posterior sampling record
  n.out <- 100
  n.record <- n.iterations - n.burnin
  out.unchanged <- floor(n.unchanged/n.out)
  
  ###Control the number of output from Cpp
  times<-ceiling(n.iterations/n.out)#times of output 
  n.last <- n.iterations - n.out * (times-1)
  
  args=list(Reads=as.integer(Read), K = as.integer(K),  
            BNum = as.integer(B), CNum = as.integer(nb), 
            GNum=as.integer(G), #observed data
            alpha.init = alpha.init, beta.init = beta.init, 
            nu.init= nu.init, delta.init = delta.init, 
            gamma.init = gamma.init, phi.init = phi.init, 
            pi.init = pi.init,p.init = p.init,tau0.init = tau0.init, #initial value 
            w.init = as.integer(w.init), z.init = as.integer(z.init), 
            x.init = as.integer(x.init), l.init = as.integer(l.init),
            xi = xi, mu.a = mu.a, mu.c = mu.c, sigma.a = sigma.a, 
            sigma.c = sigma.c, sigma.d = sigma.d, sigma.z = sigma.z, 
            kappa = kappa, tau = tau, a.tau = a.tau, b.tau = b.tau, 
            a.p = a.p, b.p = b.p, tau1 = tau1,#prior
            n.out= as.integer(n.out),ind_update_ptau = as.integer(0)) #MCMC iteration setting
  
  t.start <- Sys.time()
  
  
  for(i in seq_len(times-1)){
    results<-.Call('BZINBBUS', args)
    if(showIteration){
      message(paste0("Finish the ",i * n.out,"-th iterations",sep=" "))
    }
    
    #load the posterior sampling after n.out times iterations
    alpha_temp <- matrix(as.matrix(results[[1]]),ncol = n.out)#each row is a sample of alpha_1, alpha_2, ..., alpha_G
    beta_temp <- matrix(as.matrix(results[[2]]),ncol=n.out)#each row is a sample of beta_11,beta_21,...,beta_G1,beta_12,beta_22,...,beta_GK
    nu_temp <- matrix(as.matrix(results[[3]]),ncol=n.out)#each row is a sample of nu_11,beta_21,...,beta_B1,nu_12,nu_22,...,nu_BG
    delta_temp <- matrix(as.matrix(results[[4]]),ncol = n.out) # N * n.out
    gamma_temp <- matrix(as.matrix(results[[5]]),ncol = n.out) # B * 2 * n.out
    phi_temp <- matrix(as.matrix(results[[6]]),ncol = n.out)   # B * G * n.out
    pi_temp <- matrix(as.matrix(results[[7]]),ncol = n.out)    # K * B * n.out
    w_temp <- matrix(as.matrix(results[[8]]),ncol = n.out)     # N * n.out
    p_temp <- matrix(as.matrix(results[[9]]),ncol = n.out)     # n.out
    tau0_temp <- matrix(as.matrix(results[[10]]),ncol = n.out) # n.out
    L_temp <- matrix(as.matrix(results[[11]]),ncol = n.out)    # G * K * n.out
    
    #modify the value of non-pointers in the args
    args$p.init <- p_temp[n.out]
    args$tau0.init <- tau0_temp[n.out]
    if(i == out.unchanged){
     args$ind_update_ptau = as.integer(1) 
    }
    
    #save the posterior sampling into a txt file after burnin
    #if(n.burnin < i * n.out){#start to save
      
      # if n.burnin > (i-1) * n.out, burnin is in the inner of this output, start at n.burnin - (i-1) * n.out + 1
      # otherwise, start at 1
      #save.start <- max(1,n.burnin - (i-1) * n.out + 1)
    save.start <- 1
    write.table(t(alpha_temp[,save.start:n.out]),alpha.file,append = T,row.names = F,col.names = F)
    write.table(t(beta_temp[,save.start:n.out]),beta.file,append = T,row.names = F,col.names = F)
    write.table(t(nu_temp[,save.start:n.out]),nu.file,append = T,row.names = F,col.names = F)
    write.table(t(delta_temp[,save.start:n.out]),delta.file,append = T,row.names = F,col.names = F)
    write.table(t(gamma_temp[,save.start:n.out]),gamma.file,append = T,row.names = F,col.names = F)
    write.table(t(phi_temp[,save.start:n.out]),phi.file,append = T,row.names = F,col.names = F)
    write.table(t(pi_temp[,save.start:n.out]),pi.file,append = T,row.names = F,col.names = F)
    write.table(t(w_temp[,save.start:n.out]),w.file,append = T,row.names = F,col.names = F)
    write.table(p_temp[,save.start:n.out],prob.file,append = T,row.names = F,col.names = F)
    write.table(tau0_temp[,save.start:n.out],tau0.file,append = T,row.names = F,col.names = F)
    write.table(t(L_temp[,save.start:n.out]),l.file,append = T,row.names = F,col.names = F)
    #}
  }
  
  #the last output
  results<-.Call('BZINBBUS', args)
  if(showIteration){
    message(paste0("Finish the ", n.iterations,"-th iterations",sep=" "))
  }
  
  #load the posterior sampling after n.out times iterations
  alpha_temp <- matrix(as.matrix(results[[1]]),ncol = n.out)#each row is a sample of alpha_1, alpha_2, ..., alpha_G
  beta_temp <- matrix(as.matrix(results[[2]]),ncol=n.out)#each row is a sample of beta_11,beta_21,...,beta_G1,beta_12,beta_22,...,beta_GK
  nu_temp <- matrix(as.matrix(results[[3]]),ncol=n.out)#each row is a sample of nu_11,beta_21,...,beta_B1,nu_12,nu_22,...,nu_BG
  delta_temp <- matrix(as.matrix(results[[4]]),ncol = n.out) # N * n.out
  gamma_temp <- matrix(as.matrix(results[[5]]),ncol = n.out) # B * 2 * n.out
  phi_temp <- matrix(as.matrix(results[[6]]),ncol = n.out)   # B * G * n.out
  pi_temp <- matrix(as.matrix(results[[7]]),ncol = n.out)    # K * B * n.out
  w_temp <- matrix(as.matrix(results[[8]]),ncol = n.out)     # N * n.out
  p_temp <- matrix(as.matrix(results[[9]]),ncol = n.out)     # n.out
  tau0_temp <- matrix(as.matrix(results[[10]]),ncol = n.out) # n.out
  L_temp <- matrix(as.matrix(results[[11]]),ncol = n.out)    # G * K * n.out
  
  #modify the value of non-pointers in the args
  args$p.init <- p_temp[n.out]
  args$tau0.init <- tau0_temp[n.out]
  
  #save the posterior sampling into a txt file after burnin
  save.start <- 1
  write.table(t(alpha_temp[,save.start:n.out]),alpha.file,append = T,row.names = F,col.names = F)
  write.table(t(beta_temp[,save.start:n.out]),beta.file,append = T,row.names = F,col.names = F)
  write.table(t(nu_temp[,save.start:n.out]),nu.file,append = T,row.names = F,col.names = F)
  write.table(t(delta_temp[,save.start:n.out]),delta.file,append = T,row.names = F,col.names = F)
  write.table(t(gamma_temp[,save.start:n.out]),gamma.file,append = T,row.names = F,col.names = F)
  write.table(t(phi_temp[,save.start:n.out]),phi.file,append = T,row.names = F,col.names = F)
  write.table(t(pi_temp[,save.start:n.out]),pi.file,append = T,row.names = F,col.names = F)
  write.table(t(w_temp[,save.start:n.out]),w.file,append = T,row.names = F,col.names = F)
  write.table(p_temp[,save.start:n.out],prob.file,append = T,row.names = F,col.names = F)
  write.table(tau0_temp[,save.start:n.out],tau0.file,append = T,row.names = F,col.names = F)
  write.table(t(L_temp[,save.start:n.out]),l.file,append = T,row.names = F,col.names = F)
  
  x_temp <- matrix(args$x.init, N, G)
  write.table(x_temp,x.file,row.names = F,col.names = F)
  
  
  t.end <- Sys.time()
  
  message(paste0("  The MCMC sampling takes: ", round(difftime(t.end, t.start,
                                                                units = "mins"), 3), " mins", "\n"))
  #######################
  # Posterior inference #
  #######################
  #load posterior sampling from the "MCMC_sampling" folder
  setwd(working_dir)
  sampling_dir <- paste0(working_dir,"/MCMC_sampling",sep="")
  
  message("  loading the posterior samples...\n")
  #load(paste0(sampling_dir,"/raw_data.Rdata",sep=""))
  #load posterior sampling of cell type related parameters 
  
  t.start<-Sys.time()
  message("  calculating posterior means and posterior modes...\n")
  
  pi.file <- paste0(sampling_dir,"/pi.txt",sep="")
  pi.post<-as.matrix(read.big.matrix(pi.file,skip = n.burnin, sep=" ",type="double"))
  w.file <- paste0(sampling_dir,"/w.txt",sep="")
  w.post<-as.matrix(read.big.matrix(w.file,skip = n.burnin,sep=" ",type="integer"))
  gamma.file <- paste0(sampling_dir,"/gamma.txt",sep="")
  gamma.post<-as.matrix(read.big.matrix(gamma.file,skip = n.burnin,sep=" ",type="double"))
  alpha.file <- paste0(sampling_dir,"/alpha.txt",sep="")
  alpha.post<-as.matrix(read.big.matrix(alpha.file,skip = n.burnin,sep=" ",type="double"))
  beta.file <- paste0(sampling_dir,"/beta.txt",sep="")
  beta.post<-as.matrix(read.big.matrix(beta.file,skip = n.burnin,sep=" ",type="double"))
  delta.file <- paste0(sampling_dir,"/delta.txt",sep="")
  delta.post<-as.matrix(read.big.matrix(delta.file,skip = n.burnin,sep=" ",type="double"))
  nu.file <- paste0(sampling_dir,"/nu.txt",sep="")
  nu.post<-as.matrix(read.big.matrix(nu.file,skip = n.burnin,sep=" ",type="double"))
  phi.file <- paste0(sampling_dir,"/phi.txt",sep="")
  phi.post<-as.matrix(read.big.matrix(phi.file,skip = n.burnin,sep=" ",type="double"))
  l.file <- paste0(sampling_dir,"/l.txt",sep="")
  l.post<-as.matrix(read.big.matrix(l.file,skip = n.burnin,sep=" ",type="integer"))
  tau0.file <- paste0(sampling_dir,"/tau0.txt",sep="")
  tau0.post<-as.matrix(read.big.matrix(tau0.file,skip = n.burnin,sep=" ",type="double"))
  p.file <- paste0(sampling_dir,"/p.txt",sep="")
  p.post<-as.matrix(read.big.matrix(p.file,skip = n.burnin,sep=" ",type="double"))
  
  x.file <- paste0(sampling_dir,"/x.txt",sep="")
  x.post <- as.matrix(read.big.matrix(x.file,sep=" ",type="integer"))
  
  alpha.est<-apply(alpha.post,2,mean)
  beta.est<-matrix(apply(beta.post,2,mean),nrow=G)
  nu.est<-matrix(apply(nu.post,2,mean),nrow=B)
  delta.est<-apply(delta.post,2,mean)
  gamma.est<-matrix(apply(gamma.post,2,mean),nrow=B)
  logmu.est<-beta.est+alpha.est
  phi.est<-matrix(apply(phi.post,2,mean),nrow=B)
  pi.est<-matrix(apply(pi.post,2,mean),nrow=K)
  tau0.est<-mean(tau0.post)
  p.est <- mean(p.post)
  
  w.est <- rep(0,N)
  for(i in 1:N){
    w.est[i] <- .getmode(w.post[,i]) + 1
  }
  
  PPI.est <- apply(l.post,2,mean)
  PPI.est <- matrix(PPI.est,G,K)
  
  alpha.sd<-apply(alpha.post,2,sd)
  beta.sd<-matrix(apply(beta.post,2,sd),nrow=G)
  nu.sd<-matrix(apply(nu.post,2,sd),nrow=B)
  delta.sd<-apply(delta.post,2,sd)
  gamma.sd<-matrix(apply(gamma.post,2,sd),nrow=B)
  phi.sd<-matrix(apply(phi.post,2,sd),nrow=B)
  pi.sd<-matrix(apply(pi.post,2,sd),nrow=K)
  p.sd <- sd(p.post)
  tau0.sd <- sd(tau0.post)
  
  t.end<-Sys.time()
  message(paste0("  calculating posterior means and posterior takes: ", round(difftime(t.end, t.start,
                                                                                       units = "mins"), 3), " mins", "\n"))
  ###Generate the output "BUSseqfits" object
  #transfer Read, x.post, delta.est, w.est, delta.sd as a list
  Read_list <- list()
  Read_sim_list <- list()
  delta.est_list <- list()
  delta.sd_list <- list()
  w_list <- list()
  cell_index <- 0
  for(b in 1:B){
    Read_list[[b]] <- t(Read[cell_index + 1:nb[b],])
    Read_sim_list[[b]] <- t(x.post[cell_index + 1:nb[b],])
    delta.est_list[[b]] <- delta.est[cell_index + 1:nb[b]]
    delta.sd_list[[b]] <- delta.sd[cell_index + 1:nb[b]]
    w_list[[b]] <- w.est[cell_index + 1:nb[b]]
    cell_index <- cell_index + nb[b]
  }
  
  output <- list(CountData_raw = Read_list, CountData_imputed = Read_sim_list,
                 #dimensions
                 n.cell = N, n.gene = G, n.batch = B, 
                 n.perbatch = nb, n.celltype = K, n.record = n.record,
                 #posterior mean or mode of parameters
                 gamma.est = gamma.est, alpha.est = alpha.est, 
                 beta.est = beta.est, nu.est = nu.est,
                 delta.est = delta.est_list, phi.est = phi.est, pi.est = pi.est,
                 w.est = w_list, p.est = p.est, tau0.est = tau0.est,
                 PPI.est = PPI.est,
                 #posterior sd of pararmeters
                 alpha.sd = alpha.sd, beta.sd = beta.sd, nu.sd = nu.sd,
                 delta.sd = delta.sd_list, phi.sd = phi.sd, pi.sd = pi.sd,
                 p.sd = p.sd, tau0.sd = tau0.sd)
  
  class(output) <- "BUSseqfits"
  
  t.start<-Sys.time()
  message("  calculating BIC...\n")
  output$BIC <- .cal_BIC_BUSseq(output)
  t.end<-Sys.time()
  message(paste0("  calculating BIC takes: ", round(difftime(t.end, t.start,
                                                                            units = "mins"), 3), " mins", "\n"))
  # t.start<-Sys.time()
  # message("  correcting read counts...\n")
  # output$CountData_corrected <- adjusted_values_BUSseq(output)
  # t.end<-Sys.time()
  # message(paste0("  Correcting read counts takes: ", round(difftime(t.end, t.start,
  #                                                                  units = "mins"), 3), " mins", "\n"))
  
  return(output)
  
}

##################################
# Useful Outputs from BUSseqfits #
##################################
#obtain the cell type indicators for samples
celltypes <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .B <- BUSseqfits_obj$n.batch
    .w <- BUSseqfits_obj$w.est
    
    for(b in seq_len(.B)){
      
      message(paste0("Batch ", b, " cells' cell type indicators: ",
                .w[[b]][1],",",.w[[b]][2],",",.w[[b]][3], 
                "... ...\n"))
    }
    message("The output format is a list with length",
            " equal to the batch number.\n")
    message(paste0("Each element of the list is a cell type indicator vector in",
                   " that batch.\n"))
    return(.w)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the dropout intercept and odds ratio
dropout_coefficient_values <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .gamma<-BUSseqfits_obj$gamma.est
    message("The output format is a matrix.\n")
    message(paste0("Each row represents a batch, the first column corresponds",
                   " to intercept and the second column is the odd ratio.\n"))
    return(.gamma)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the log-scale baseline expression values
baseline_expression_values <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .alpha<-BUSseqfits_obj$alpha.est
    message("The output format is a vector.\n")
    return(.alpha)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the cell type effects
celltype_effects <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .beta<-BUSseqfits_obj$beta.est
    message("The output format is a matrix.\n")
    message(paste0("Each row represents a gene, and each column corresponds",
                   " to a cell type.\n"))
    return(.beta)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the cell-tpye-specific mean expression levels
celltype_mean_expression <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .alpha<-BUSseqfits_obj$alpha.est
    .beta<-BUSseqfits_obj$beta.est
    
    mu <- exp(.alpha+.beta)
    
    message("The output format is a matrix.\n")
    message(paste0("Each row represents a gene, and each column corresponds",
                   " to a cell type.\n"))
    return(mu)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the location batch effects
location_batch_effects <- function(BUSseqfits_obj){
  
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .nu<-BUSseqfits_obj$nu.est
    message("The output format is a matrix.\n")
    message("Each row represents a batch, and each column",
            " corresponds to a gene.\n")
    return(.nu)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the scale batch effects
overdispersions <- function(BUSseqfits_obj){
  
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .phi<-BUSseqfits_obj$phi.est
    message("The output format is a matrix.\n")
    message("Each row represents a batch, and each column",
            " corresponds to a gene.\n")
    return(.phi)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the cell-specific size effects
cell_effect_values <- function(BUSseqfits_obj){
  
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .delta<-BUSseqfits_obj$delta.est
    message("The output format is a list with length equal to the batch number.\n")
    message(paste0("Each element of the list is a cell-specific size factor vector of",
                   " that batch.\n"))
    return(.delta)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}


#obtain the instrinsic genes
intrinsic_genes_BUSseq <- function(BUSseqfits_obj,fdr_threshold=0.05){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    
    .PPI<-BUSseqfits_obj$PPI.est
    postprob_DE_threshold <- .postprob_DE_thr_fun(.PPI,fdr_threshold)
    postprob_DE_threshold <- max(postprob_DE_threshold, 0.5)
    L.est <- .estimate_IG_indicators(.PPI,postprob_DE_threshold)
    D.est <- .IG_index(L.est)
    return(D.est)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the BIC score
BIC_BUSseq <- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    .BIC<-BUSseqfits_obj$BIC
    message("BIC is ", .BIC, "\n")
    message("The output is a scalar.\n")
    return(.BIC)
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the raw count data
raw_read_counts<- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    
    Read_raw <- BUSseqfits_obj$CountData_raw
    
    class(Read_raw) <- "CountData"
    
    message("The output format is a \"CountData\" object with length equal to the batch number.\n")
    message("Each element of the object is the raw read count matrix.\n")
    message("In each matrix, each row represents a gene and each column correspods to a cell.\n")
    
    return(Read_raw)
    
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obtain the underlying true count data
imputed_read_counts<- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    
    Read_imputed <- BUSseqfits_obj$CountData_imputed
    
    class(Read_imputed) <- "CountData"
    
    message("The output format is a \"CountData\" object with length equal to the batch number.\n")
    message("Each element of the object is the imputed read count matrix.\n")
    message("In each matrix, each row represents a gene and each column correspods to a cell.\n")
    
    return(Read_imputed)
    
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

#obatin the corrected count data
corrected_read_counts<- function(BUSseqfits_obj){
  if(is(BUSseqfits_obj,"BUSseqfits")){
    
    Truereads_list <- BUSseqfits_obj$CountData_imputed
    .B <- BUSseqfits_obj$n.batch
    .nb <- BUSseqfits_obj$n.perbatch
    .N <- BUSseqfits_obj$n.cell
    .G <- BUSseqfits_obj$n.gene
    .K <- BUSseqfits_obj$n.celltype
    .gamma <- BUSseqfits_obj$gamma.est
    .logmu <- BUSseqfits_obj$alpha.est + BUSseqfits_obj$beta.est
    .nu <- BUSseqfits_obj$nu.est
    .delta <- BUSseqfits_obj$delta.est
    .phi <- BUSseqfits_obj$phi.est
    .w <- BUSseqfits_obj$w.est
    
    Truereads <- NULL
    delta_vec <- NULL
    w_vec <- NULL
    for(b in 1:.B){
      Truereads <- rbind(Truereads,t(Truereads_list[[b]]))
      delta_vec <- c(delta_vec, .delta[[b]])
      w_vec <- c(w_vec, .w[[b]])
    }
    
    t.start<-Sys.time()
    message("  correcting read counts...\n")
    
    read_corrected <- Truereads
    for(b in 1:.B){
      for(i in 1:.nb[b]){
        if(b==1){
          row.index <- i
        }else{
          row.index <- i + sum(.nb[1:(b-1)])
        }
        for(j in 1:.G){
          p_x<-pnbinom(Truereads[row.index,j],size = .phi[b,j], mu = exp(.logmu[j,w_vec[row.index]] + .nu[b,j] + delta_vec[row.index]))
          p_xminus1<-pnbinom(Truereads[row.index,j]-1,size = .phi[b,j], mu = exp(.logmu[j,w_vec[row.index]] + .nu[b,j] + delta_vec[row.index]))
          u <- runif(1,min=p_xminus1,max=p_x)
          u <- min(u,0.999)
          read_corrected[row.index,j] <- qnbinom(u,size =.phi[1,j], mu = exp(.logmu[j,w_vec[row.index]]))
        }
        #print(paste("Finish the ",row.index,"-th cell",sep=""))
      }
    }
    
    Read_corrected <- list()
    cell_index <- 0
    for(b in 1:.B){
      Read_corrected[[b]] <- t(read_corrected[cell_index + 1:.nb[b],])
      cell_index <- cell_index + .nb[b]
    }
    
    class(Read_corrected) <- "CountData"
    
    t.end<-Sys.time()
    message(paste0("  Correcting read counts takes: ", round(difftime(t.end, t.start,
                                                                      units = "mins"), 3), " mins", "\n"))
    
    message("The output format is a \"CountData\" object with length equal to the batch number.\n")
    message("Each element of the object is the corrected read count matrix.\n")
    message("In each matrix, each row represents a gene and each column correspods to a cell.\n")
    
    return(Read_corrected)
    
  }else{
    stop("BUSseqfits_obj must be  a \"BUSseqfits\" object!\n")
  }
}

##############################################################################
# print and summary
##############################################################################
#print BUSseqfits
print.BUSseqfits <- function(x, ...){
  BUSseqfits <- x
  .G <- BUSseqfits$n.gene
  .B <- BUSseqfits$n.batch
  
  
  message("Cell type indicators:\n")
  .w <- BUSseqfits$w.est
  .nb <- BUSseqfits$n.perbatch
  for(b in seq_len(.B)){
    
    message(paste0("Batch ", b, " cells' cell type indicators: ",
                   .w[[b]][1],",",.w[[b]][2],",",.w[[b]][3], 
                   "... ...\n"))
  }
  message("\n")
  
  message("The estimated location batch effects:\n")
  .nu <- BUSseqfits$nu.est
  for(b in seq_len(.B)){
    message(paste0("   Batch ", b, " location batch effects are: ",
          .nu[b,1],",",.nu[b,2],",",.nu[b,3], "... ...\n"))
  }
  message("\n")
  
  message("The estimated overdispersions:\n")
  .phi <- BUSseqfits$phi.est
  for(b in seq_len(.B)){
    message(paste0("   Batch ", b, " scale batch effects are: ",
          .phi[b,1],",",.phi[b,2],",",.phi[b,3], "... ...\n"))
  }
  message("\n")
}

#summarize BUSseqfits
summary.BUSseqfits <- function(object, ...){
  BUSseqfits <- object
  .G <- BUSseqfits$n.gene
  .B <- BUSseqfits$n.batch
  .K <- BUSseqfits$n.celltype
  .N <- BUSseqfits$n.cell
  
  num_records <- BUSseqfits$n.record
  message(c("B = ", .B, " batches\n"))
  message(c("G = ", .G, " genes\n"))
  message(c("K = ", .K, " cell types\n"))
  message(c("N = ", .N, " cells in total\n"))
  message(c("n.record = ", num_records," iterations are recorded.\n\n"))
  message("BUSseqfits is an R list that contains the following main elements:\n\n")
  message(paste0("   BUSseqfits$w.est : the estimated cell type indicators,",
             " a list with length equal to B.\n"))
  message(paste0("   BUSseqfits$pi.est : the estimated cell type proportions across batches,",
             " a K by B matrix.\n"))
  message(paste0("   BUSseqfits$gamma.est : the estimated the coefficients of the logistic ",
                 "regression for the dropout events, a B by 2 matrix\n"))
  message(paste0("   BUSseqfits$alpha.est : the estimated log-scale baseline expression levels,",
             " a vector with length G.\n"))
  message(paste0("   BUSseqfits$beta.est : the estimated cell type effects,",
             " a G by K matrix.\n"))
  message(paste0("   BUSseqfits$delta.est : the estimated cell-specific effects,",
                 " a list with length equal to B.\n"))
  message(paste0("   BUSseqfits$nu.est : the estimated location batch effects,",
                 " a B by G matrix.\n"))
  message(paste0("   BUSseqfits$phi.est : the estimated overdispersion parameters,",
             " a B by G matrix.\n"))
  message("   BUSseqfits$BIC : the BIC, a scalar.\n")
  message(paste0("   BUSseqfits$PPI.est : the posterior marginal probablity of being ",
             "an intrinsic gene, a G by K matrix.\n"))
  message("   For more output values, please use \"?BUSseq_MCMC\"\n")
  message("\n")
}

#print CountData
print.CountData <- function(x, ...){
  CountData <- x

  .B <- length(CountData)
  .G <- nrow(CountData[[1]])
  
  message(paste0("There are ", .B, " batches and ", .G, " genes.\n"))
  for(b in seq_len(.B)){
    .nperbatch <- ncol(CountData[[b]])
    message(paste0("Batch ", b, " contains ", .nperbatch, " cells, and their read counts in all genes are: \n"))
    message(paste0("Gene 1: ", CountData[[b]][1,1],", ", CountData[[b]][1,2], ", ", CountData[[b]][1,3], ", ... ...\n"))
    message(paste0("Gene 2: ", CountData[[b]][2,1],", ", CountData[[b]][2,2], ", ", CountData[[b]][2,3], ", ... ...\n"))
    message(paste0("Gene 3: ", CountData[[b]][3,1],", ", CountData[[b]][3,2], ", ", CountData[[b]][3,3], ", ... ...\n"))
    message("   ... ...\n\n")
  }
  message("\n")
  
}

#summarize CountData
summary.CountData <- function(object, ...){

  CountData <- object
  
  .B <- length(CountData)
  .G <- nrow(CountData[[1]])
  
  message(paste0("There are ", .B, " batches and ", .G, " genes.\n"))
  for(b in seq_len(.B)){
    .nperbatch <- ncol(CountData[[b]])
    message(paste0("Batch ", b, " contains ", .nperbatch, " cells."))
  }
  message("\n")
}

########################################################################
# Visualization
########################################################################
#visualize the read counts data by stacking all gene expression matrices 
heatmap_data_BUSseq <- function(CountData_obj, gene_set = NULL, project_name="BUSseq_heatmap", 
                                    image_dir = NULL, color_key_seq = NULL, 
                                    image_width = 1440, image_height = 1080){
  
  if(is(CountData_obj,"CountData")){
  .B <- length(CountData_obj)
  .G <- nrow(CountData_obj[[1]])
  .nb <- rep(NA,.B)
  for(b in seq_len(.B)){
    .nb[b] <- ncol(CountData_obj[[b]])
  }

  if(is.null(gene_set)){
    gene_set <- 1:.G
  }

  #heatmap cell colors
  colfunc <- colorRampPalette(c("grey", "black"))
  #batch colors
  color_batch_func <- colorRampPalette(c("#EB4334","#FBBD06","#35AA53","#4586F3"))
  
  color_batch <- color_batch_func(.B)
  
  color_batch2 <- NULL
  
  for(b in seq_len(.B)) {
    color_batch2 <- c(color_batch2, rep(color_batch[b], .nb[b]))
  }
  log1p_mat <- NULL
  for(b in seq_len(.B)){
    log1p_mat <- cbind(log1p_mat, log1p(CountData_obj[[b]]))
  }
  log1p_mat_interest <- log1p_mat[gene_set, ]
  
  if(is.null(color_key_seq)){
    range_data <- range(log1p_mat_interest)
    color_key_seq <- seq(from = floor(range_data[1]) - 0.5, to = ceiling(range_data[2]) + 0.5, length.out = 11)
  }
  
  if(is.null(image_dir)){
    image_dir <- "./image"
  }
  #create the folder
  dir.create(image_dir,showWarnings = F)
  
  png(paste(image_dir,"/",project_name,"_log1p_data.png",sep=""),width = image_width, height = image_height)
  heatmap.2(log1p_mat_interest,
            dendrogram = "none",#with cluster tree
            Rowv = FALSE, Colv = FALSE,
            #xlab = "SUBTYPE", ylab = "GENE", #main = "log(True subtype mean)",
            labRow = FALSE, labCol = FALSE,
            ColSideColors = color_batch2,
            #RowSideColors = genetype_color,
            col = colfunc(length(color_key_seq)-1),breaks = color_key_seq,
            density.info="histogram",
            hclustfun = function(c)hclust(c,method="average"),keysize = 0.8, cexRow=0.5,trace = "none")#font size
  dev.off()
  }else{
    stop("CountData_obj must be  a \"CountData\" object!\n")
  }
}
