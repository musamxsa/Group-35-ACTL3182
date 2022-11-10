library("quantmod")
library("tidyquant")
library("imputeTS")
library("dplyr")
library("rstan")
library("tidyr")
library(tibble)
library("readr")
library("quadprog")





Tickers=c("^GSPC", "^DJI","^GSPE","^SP500-35","^SP500-45","^SP500-60",
          "SHY","^SPGSCI" )



getSymbols(Tickers, from ='2002-07-30', to ='2022-11-02')

OG_GSPC<-c(GSPC$GSPC.Adjusted)
OG_DJI<-c(DJI$DJI.Adjusted)
OG_GSPE<-c(GSPE$GSPE.Adjusted)
OG_SP35<-c(`SP500-35`$`SP500-35.Adjusted`)
OG_SP45<-c(`SP500-45`$`SP500-45.Adjusted`)
OG_SP60<-c(`SP500-60`$`SP500-60.Adjusted`)
OG_SHY<-c(SHY$SHY.Adjusted)
OG_SPGSCI<-c(SPGSCI$SPGSCI.Adjusted)


###############

#For log returns

Ret_Log<- function(x){

seq<-c()  
    
  for (i in 2:length(x)){
    
    d<- log(as.numeric(x[i])/as.numeric(x[i-1]))
    seq<-c(seq,d)
    
    seq[is.na(seq)]<-0
    
  }
  
  return(seq)
  
}


ret<-data.frame(
  
  LOG_GSPC<-Ret_Log(OG_GSPC),
  LOG_DJI<-Ret_Log(OG_DJI),
  LOG_GSPE<-Ret_Log(OG_GSPE),
  LOG_SP35<-Ret_Log(OG_SP35),
  LOG_SP45<-Ret_Log(OG_SP45),
  LOG_SP60<-Ret_Log(OG_SP60),
  LOG_SHY<-Ret_Log(OG_SHY),
  LOG_SPGSCI<-Ret_Log(OG_SPGSCI) 
  
  
)


mean_ret<-apply(ret, 2, mean)

mean_ret


cov_mat <- cov(ret)

cov_mat

##### 
T <- nrow(ret)
N <- ncol(ret)
nu <- 12
data_stan <- list(
  T = T,
  N = N,
  nu = nu,
  tau = 200,
  eta = mean_ret,
  R = as.matrix(ret),
  omega = cov_mat * (nu - N -1)
)

fit <- stan(
  file = "bay_port.stan",
  data = data_stan,
  chains = 4,
  warmup = 300,
  iter = 500,
  cores = 2
)

saveRDS(fit, 'stan_fit.rds')
fit <-readRDS('stan_fit.rds')

traceplot(fit, nrow = 4, pars = c('mu', 'sigma'))
print(fit, pars = c('mu', 'sigma'))

#Extract draws from the posterior
list_of_draws <- extract(fit)
#draws from the posterior distribution of sigma
sigma_post <- list_of_draws$sigma
#draws from the posterior distribution of mu
mu_post <- list_of_draws$mu

fit <-readRDS('stan_fit.rds')

#Some diagnostics
traceplot(fit, nrow = 4, pars = c('mu', 'sigma'))

#Extract draws from the posterior
list_of_draws <- extract(fit)
#draws from the posterior distribution of sigma
sigma_post <- list_of_draws$sigma
#draws from the posterior distribution of mu
mu_post <- list_of_draws$mu

#Solves the optimization
#problem for distinct
#draws of the posterior
#parameters.
#Each draw will have an
#efficient frontier associated
#with it
n_frontiers <- 10
set.seed(54321)
n_draws <- nrow(mu_post)
sample_index <- sample(n_draws,
                       size = n_frontiers,
                       replace = FALSE)

#target (annual) return
mu_target <- seq(0.02,
                 max(apply(mu_post, 2, max)),
                 le = 100)

#Number of assets
N <- ncol(mu_post)
N
#aux for detecting last
#frontier
aux_last <- 0

par(bg = '#EEEEEC',
    mfcol = c(ceiling(n_frontiers / 2),2))
for(idx in sample_index){
  aux_last <- aux_last + 1
  
  #for the table that will
  #be use to create the graph
  mu_tib <-  c()
  var_opt <- c()
  #Last frontier is classical
  #framework
  if(aux_last == 1){
    mu_draw <- mean_ret
    sig_draw <- cov_mat
  }
  else{
    #draw from mu
    mu_draw <- mu_post[idx,]
    
    #draw from sigma
    sig_draw <- sigma_post[idx, ,]
  }
  
  #Solves the optimization
  #problem for each target value
  for(val in mu_target){
    
    A <- matrix(0, nrow = N,ncol = 2)
    #sum of weights equals 1
    A[,1] <- 1
    
    #the target return constrain
    A[,2] <- mu_draw * 252
    
    b0 <- c(1, val)
    sol <- solve.QP(2 * 252 * sig_draw,
                    dvec = rep(0, N),
                    Amat = A,
                    bvec = b0,
                    meq = 2)
    mu_tib <- c(mu_tib, val)
    var_opt <- c(var_opt, sol$value )
  }
  
  #Creates the tibble
  #For the chart
  chart_tib <- tibble(mu =  100 * mu_tib,
                      std = 100 * sqrt(var_opt))
  
  main <- 'Efficient Frontier'
  
  if(aux_last == 1){
    main <- 'Efficient Frontier \n Classical Framework'
  }
  
  plot(chart_tib$std, chart_tib$mu,
       col = 'blue', lwd = 2.5,
       main = main,
       sub = 'Annualized figures',
       xlab = 'Standard Deviation %',
       ylab = 'Expected return %',
       type = 'l',
       ylim = c(100 * min(mu_tib), 100 * max(mu_tib)))
  grid(col = 'black', lwd = 1.5)
}


#Solves the optimization
#and averages the solutions
#found for every target level
#n_frontiers <- nrow(mu_post) # the whole data
n_frontiers <- 200


max(apply(mu_post,2,max))

#target (annual) return
mu_target <- seq(0.02,
                 1.118302,
                 le = 100)

#This tibble will store the
#optimal weights
opt_weights <- tibble(GSPC = NULL,
                      DJI = NULL,
                      GSPE = NULL,
                      SP35 = NULL,
                      SP45 = NULL,
                      SP60 = NULL,
                      SHY = NULL,
                      SPGSCI = NULL,
                      target = NULL,
                      sd = NULL) 

#Number of assets
N <- ncol(mu_post)
N
for(idx in 1:n_frontiers){
  
  #for the table that will
  #be use to create the graph
  mu_tib <-  c()
  var_opt <- c()
  
  #draw from mu
  mu_draw <- mu_post[idx,]
  
  #draw from sigma
  sig_draw <- sigma_post[idx, ,]
  
  #Solves the optimization
  #problem for each target value
  for(val in mu_target){
    
    A <- matrix(0, nrow = N,ncol = 2)
    #sum of weights equals 1
    A[,1] <- 1
    
    #the target return constrain
    A[,2] <- mu_draw * 252
    
    b0 <- c(1, val)
    sol <- solve.QP(2 * 252 * sig_draw,
                    dvec = rep(0, N),
                    Amat = A,
                    bvec = b0,
                    meq = 2)
    var_opt <- c(var_opt, sol$value )
    
    #updates opt_weigths tibble
    row <- c(sol$solution, val, sol$value)
    opt_weights <- rbind(opt_weights, row)
  }
  
}
names(opt_weights) <- c('GSPC',
                        'DJI',
                        'GSPE',
                        'SP35',
                        'SP45',
                        'SP60',
                        'SHY',
                        'SPGSCI',
                        'target',
                        'std')

write_csv(opt_weights, "opt_weights.csv") 


#tibble with the optimal
#weights
opt_weights <- read_csv('opt_weights.csv')

#Averages the values for each
#target value
mean_opt_w <- opt_weights %>%
  group_by(target) %>%
  mutate(mean_GSPC = mean(GSPC),
         mean_DJI = mean(DJI),
         mean_GSPE = mean(GSPE),
         mean_SP35 = mean(SP35),
         mean_SP45 = mean(SP45),
         mean_SP60 = mean(SP60),
         mean_SHY = mean(SHY),
         mean_SPGSCI = mean(SPGSCI),
         mean_std = mean(std)) %>%
  select(mean_GSPC,
         mean_GSPE,
         mean_DJI,
         mean_SP35,
         mean_SP45,
         mean_SP60,
         mean_SHY,
         mean_SPGSCI,
         target,
         mean_std) %>%
  unique()

par(bg = '#EEEEEC',
    mfcol = c(1,1))

main <- 'Average Efficient Frontier'
plot(100 * mean_opt_w$mean_std,
     100 * mean_opt_w$target,
     col = 'blue', lwd = 2.5,
     main = main,
     sub = 'Annualized figures',
     xlab = 'Standard Deviation %',
     ylab = 'Expected return %',
     type = 'l')
grid(col = 'black', lwd = 1.5)

mean_opt_w$mean_GSPE






############################################################################

long <- mean_opt_w %>%
  filter((mean_GSPC >= 0 & mean_GSPC <= 1)&
  (mean_GSPE >= 0 & mean_GSPE <= 1)&
  (mean_SP35>=0 & mean_SP35<=1)&
  (mean_SP45>=0 & mean_SP45<=1)&
  (mean_SP60>=0 & mean_SP60<=1)&
  (mean_SHY>=0 & mean_SHY<=1) &
  (mean_SPGSCI>=0 & mean_SPGSCI<=1)&
  (mean_DJI>=0 & mean_DJI<=1))
  

short <- mean_opt_w %>%
  filter((mean_GSPC >= -1 & mean_GSPC <= 1) &
    (mean_DJI >= -1 & mean_DJI <= 1) &
    (mean_GSPE >= -1 & mean_GSPE <= 1)&
    (mean_SP35>=-1 & mean_SP35<=1)&
    (mean_SP45>=-1 & mean_SP45<=1)&
    (mean_SP60>=-1 & mean_SP60<=1)&
    (mean_SHY>=-1 & mean_SHY<=1)&
    (mean_SPGSCI>=-1 & mean_SPGSCI<=1))
















##### Ignore this #####
short
par(bg = '#EEEEEC',
    mfcol = c(1,2))

main_long <- 'Average Efficient Frontier \n Long Only'
main_short <- 'Average Efficient Frontier \n Short-Long'

plot(100 * long$mean_std,
     100 * long$target,
     col = 'blue', lwd = 2.5,
     main = main_long,
     sub = 'Annualized figures',
     xlab = 'Standard Deviation %',
     ylab = 'Expected return %',
     type = 'l')
grid(col = 'black', lwd = 1.5)

plot(100 * short$mean_std,
     100 * short$target,
     col = 'blue', lwd = 2.5,
     main = main_short,
     sub = 'Annualized figures',
     xlab = 'Standard Deviation %',
     ylab = 'Expected return %',
     type = 'l')
grid(col = 'black', lwd = 1.5)



