library(readxl)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

set.seed(1997)

#### Global Parameters ####

Num_sim = 10
len_path = 100

S0 = 100 # Initial stock price
q = 0 # No divident payment
r = 0.02 # The interest rates

vol = 0.05

#### Black Scholes ####
for(i in 1:10){
  St_vec = c(S0)
  #All iid normal variables are simulated
  eps = rnorm(len_path)
  
  for(j in 1:len_path){
    St_vec[j+1] = St_vec[j]+St_vec[j]*(r/len_path+sqrt(vol/len_path)*eps[j])
  }
  if(i==1){
    plot(St_vec, type = "l", ylim = c(65,150), lwd = 1.5, main = "Black-Scholes")
  }else{
    lines(St_vec, type = "l", lwd = 1.5)
  }
}


#### Heston ####
kappa = 0.5
nu = 0.05
theta = 0.2
rho = -0.75

for(i in 1:10){
  St_vec = c(S0)
  vol_vec = c(vol)
  #All iid normal variables are simulated
  eps1 = rnorm(len_path)
  eps2 = rnorm(len_path)
  
  epsS = eps1
  epsv = rho*eps1 + sqrt(1-rho^2)*eps2
  
  for(j in 1:len_path){
    St_vec[j+1] = St_vec[j]+St_vec[j]*(r/len_path+sqrt(vol_vec[j]/len_path)*epsS[j])
    vol_vec[j+1] = vol_vec[j] + (kappa*( nu-vol_vec[j] ) - (theta^2)/4 )/len_path + 
      theta*sqrt(vol_vec[j]/len_path)*epsv[j]+(theta^2)*(epsv[j]^2)/(4*len_path)
  }
  if(i==1){
    plot(St_vec, type = "l", ylim = c(65,150), lwd = 1.5, main ="Heston")
  }else{
    lines(St_vec, type = "l", lwd = 1.5)
  }
}

plot(vol_vec, type = "l", lwd = 1.5, main="Vol Heston")
