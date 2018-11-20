##Toy Velocity with unknown number of components
#library(MCMCpack)
t_0 = 1
num_samples = 1000


iterations = 2000
burn_in_percent = 0.5
end_burn = iterations*burn_in_percent




flag_sim = 1

true_lambda_k = 2 #Assume the number of mixtures follows a poisson distribution with mean lambda_k
true_gamma_pi = 1

true_k = 2 #3#max(rpois(1, lambda = true_lambda_k),1)
true_alpha_vec = c(10*100,1*50)
true_beta_vec = c(100, 50)
true_pis = c(0.6,0.4)


mix_means = true_alpha_vec/true_beta_vec; print(mix_means)
mix_var = true_alpha_vec/true_beta_vec^2; print(mix_var)





#=========================================================================
# Simualting the fake velocity data
#=========================================================================
if(flag_sim ==1){

xx_list = as.list(seq(0.0001, 10, length = 200))

true_state_vec = rmultinom(n = num_samples,size =1 , prob = true_pis )
num_per_mix = apply(true_state_vec, 1, sum)


velocity_samp = c()
for(kk in 1:true_k){
	velocity_samp = c(velocity_samp, rgamma(n = num_per_mix[kk], shape = true_alpha_vec[kk], rate= true_beta_vec[kk]))	
}
par(mfrow = c(1,3), oma = c(0,0,1.5,0))
plot(unlist(xx_list),unlist(lapply(xx_list, mixture_density_i, pis = true_pis, alphas= true_alpha_vec, betas = true_beta_vec)), type = "l", main = "True Gamma Mixture Density", ylab = "Density")

hist(velocity_samp, main = "Histogram", xlab = "Velocity");plot(ecdf(velocity_samp), main = "Emperical CDF", xlab = "Velocity")
title("Emperical Velocity Distributions", outer = TRUE)

}


#====================================================================
# Initializing Functions
#====================================================================
k_n = 3 
pis_n = c(0.5, 0.3, 0.2)
alpha_vec_n =rgamma(k_n, mean(velocity_samp)^2, var(velocity_samp))
beta_vec_n = rgamma(k_n, mean(velocity_samp), var(velocity_samp))



#Hyper parameters for pis
hyper_pis_gamma =1


#Hyper parameters for alpha
hyper_alpha_a=(mean(velocity_samp)^2/var(velocity_samp))^(2/3);
hyper_alpha_b =1/sqrt(hyper_alpha_a) 
 
 
#Hyper parameters for beta
hyper_beta_a=(mean(velocity_samp)/var(velocity_samp))^(2/3);
hyper_beta_b=1/sqrt(hyper_beta_a) 

 
#Hyper parameters for k
k_max = 20
hyper_k = 1
prior_k_prob = dpois(0:k_max, lambda =hyper_k)/ sum(dpois(0:k_max, lambda =hyper_k))

##Hyper parameters for birth/death process
birth_rate = hyper_k # This way the ratio of the prior for k in computing the death rate is just equal to 1




#====================================================================
# Output
#====================================================================
k_sample = k_all =  vector()
pi_sample = beta_sample = alpha_sample = list()


#====================================================================
# Algorithm is a go
#====================================================================
accept_count = 0


for(n in 1:iterations){
#====================================================================
# Birth and death Algorithm Steps
#====================================================================

###------------------Calculating the death rate ---------
if(k_n == 1){death_rates_j = 0
	}else{
	death_rates_j = death_rate_vec(velocity_samp, pis_n, alpha_vec_n, beta_vec_n)
	}

infinite_death = as.numeric(is.infinite(death_rates_j))
if(sum(infinite_death)>0){
	death_rates_j[which(infinite_death ==1)] = gamma(171)
}


###------------------Generating the time to next jump ---------


time_jump = rexp(1, birth_rate + sum(death_rates_j))

if(time_jump < t_0){
prob_birth = birth_rate/(birth_rate + sum(death_rates_j))

if(runif(n = 1, min = 0, max = 1) < prob_birth && k_n < k_max){
	new_k = k_n + 1
	new_pi = rbeta(n = 1, shape1 = 1, shape2 = new_k-1)
	new_beta = rgamma(n = 1, shape = hyper_beta_a, rate = hyper_beta_b)
	new_alpha = rgamma(n = 1, shape = hyper_alpha_a, rate = hyper_alpha_b)
	
	alpha_vec_n = c(alpha_vec_n, new_alpha)
	beta_vec_n = c(beta_vec_n, new_beta)
	pis_n = c(pis_n*(1 - new_pi), new_pi)
	k_n = new_k
	
}else{
	
	death_component = sample(1:k_n, size = 1, prob = death_rates_j/sum(death_rates_j))
	new_k = k_n -1
	
	alpha_vec_n = alpha_vec_n[-death_component]
	beta_vec_n = beta_vec_n[-death_component]
	pis_n = pis_n[-death_component]/(1 - pis_n[death_component])
	
	k_n = new_k	
	
	} #Updating the vectors if birth death occurs


}
	
#print(paste("mcmc step, k = ", k_n))	

#====================================================================
#  MCMC STEP
#====================================================================

    z=matrix(NA,k_n,length(velocity_samp))
 for (i in 1:length(velocity_samp)){
	z[,i]= rmultinom(1,size=1,prob= (pis_n*dgamma(velocity_samp[i], alpha_vec_n, beta_vec_n))/sum(pis_n*dgamma(velocity_samp[i], alpha_vec_n,beta_vec_n)))
    }
    
    
    
   
    nn=vector()
    for (kk in 1:k_n){
        nn[kk]=sum(z[kk,])
        beta_vec_n[kk]=rgamma(1,hyper_beta_a+nn[kk]*alpha_vec_n[kk], hyper_beta_b+sum(velocity_samp*z[kk,]))
        pr=((beta_vec_n[kk]*velocity_samp)*z[kk,])
        
        alpha_cand=rgamma(1,hyper_alpha_a,hyper_alpha_b)
        p =exp(nn[kk]*(lgamma(alpha_vec_n[kk])-lgamma(alpha_cand))+(alpha_cand-alpha_vec_n[kk])*log(prod(pr[pr>0])))
        if (p>=runif(1)){alpha_vec_n[kk]=alpha_cand;acept_count = accept_count + 1} # Acceptance rate
        }
    pis_n=rdirichlet(1,1+nn)



k_all = c(k_all, k_n)
if(n >= end_burn){
	n_index = n-end_burn +1
	
	pi_sample[[n_index]] = pis_n 
	beta_sample[[n_index]] = beta_vec_n
	alpha_sample[[n_index]] = alpha_vec_n
	k_sample = c(k_sample, k_n)
} #End of save values

	
	
	
	
}	
	


par(mfrow = c(1,3))	
plot(sample(k_sample,length(k_sample)/5, replace = F),type="l",xlab="iterations",ylab="k1",ylim=c(1,max(k_sample)+1),main="")

plotk(k_sample)
fmix(velocity_samp,k_sample,pi_sample,alpha_sample,beta_sample)
samp_seq = seq(min(velocity_samp)*0.9,max(velocity_samp)*1.2,0.01)
true_mix = unlist(lapply(as.list(samp_seq),mixture_density_i, pis = true_pis, alphas = true_alpha_vec, betas = true_beta_vec))
lines(samp_seq, true_mix, col= "blue", type="l",lty=2,lw=3)




