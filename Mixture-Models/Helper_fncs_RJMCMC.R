
require(coda) 
require(MCMCpack)

#=========================================================================
# Priors Functions
#=========================================================================

## ---------------------- Prior on k --------------------------------------
## Truncated pois(hyper_k) with k_max 


# # prior_k = function(my_k, k_max, k_prob){
	# ind_function = as.numeric(my_k <= k_max)
	# return(ind_function*k_prob[my_k + 1])	
# }


## ------------------ Prior on Mixture Probabilities ---------------------
## Symmetric Dirchlet

prior_pis  = function(my_pis, hyper_gamma){
	return(ddirichlet(my_pis, hyper_gamma) )	
}

## ------------------ Prior on Velocity alpha ---------------------
## Simple case assume they are all the same and known



## ------------------ HyperPrior for Betas ---------------------
## Gamma


prior_beta = function(my_beta, hyper_beta_a, hyper_beta_b){
	return(dgamma(my_beta, hyper_beta_a, hyper_beta_b))
}



#=========================================================================
# List Functions
#=========================================================================

sum_index_list = function(mydata, index_to_sum){sum(mydata[index_to_sum])}
prod_index_list = function(mydata, index_to_prod){prod(mydata[index_to_prod])}


simulate_gamma  = function(num_sim, parameters){rgamma(n = num_sim, parameters[1], parameters[2])}


like_list_function = function(mydata, parameters){
	return(parameters[1]*dgamma(mydata, shape = parameters[2], rate = parameters[3] ))	
}

which_equals = function(my_vector, my_number){which(my_vector == my_number)}

#=========================================================================
# Likelihood Functions
#=========================================================================


mixture_density_i = function(x,pis, alphas, betas){	
	ll_x = pis*dgamma(x, alphas, betas)
	return((sum(ll_x)))	
}

mix_like  = function(mydata, pis, alphas, betas){
	return(sum(unlist(lapply(as.list(mydata), my_mixture_density, pis = pis, alphas = alphas, betas = betas)	)))
}



#----------------------------------------------------------------------------------------------

ll_mixture_density_i = function(x,pis, alphas, betas){	
	ll_x = pis*dgamma(x, alphas, betas)
	return(log(sum(ll_x)))	
}


mix_ll = function(mydata, pis, alphas, betas){
	
	to_return = sum(unlist(lapply(as.list(mydata), ll_mixture_density_i,pis = pis, alphas = alphas, betas = betas )))	
	return(to_return)	
}



#=========================================================================
# Birth and death rate Functions
#=========================================================================


death_rate_j = function(mydata,jth,pis, alphas, betas){
	
	my_k = length(pis)
	pi_tilde = pis[-jth]/(1-pis[jth])
	alpha_tilde = alphas[-jth]
	beta_tilde = betas[-jth]
	ratio_ll = mix_ll(mydata, pi_tilde, alpha_tilde, beta_tilde) - mix_ll(mydata, pis, alphas, betas)
	#ratio_prior_k = 1 # if birth rate if lambda_b = hyperprior for k
	
	return(exp(ratio_ll))
}



death_rate_vec = function(mydata,pis, alphas, betas){
return(unlist(lapply(as.list(1:length(pis)), death_rate_j,mydata = mydata, pis = pis,alphas = alphas, betas = betas)))
}


#=========================================================================
# state_vector_prob_j
#=========================================================================

state_vec_probability_ith = function(x_i, pis, alphas, betas){
	to_return = pis*dgamma(x_i, shape = alphas, rate =betas)/sum(pis*dgamma(x_i, shape = alphas, rate =betas))
	return(to_return)	
}



#=========================================================================
# Plotting the posterior samples
#=========================================================================
plotk=function(k){
     y=vector()
     for(i in 1:max(k)){
        y[i]=sum(k[k==i])/(i*length(k))
     }
plot(x=1:max(k),y,type="h",main="",ylab="Pr(k|data)",xlab="k")
}


fmix=function(x,k,pa,a,b){
     hist(x,prob=T,xlim=c(min(x)*0.9,max(x)*1.2),ylim=c(0,.6),main="",xlab="",nclass=25,col="8")
     tt=seq(min(x)*0.9,max(x)*1.2,length=500)
     f=0*tt
     for(i in 1:length(tt)){
         for(j in 1:length(k)){
             f[i]=f[i]+sum(pa[[j]]*dgamma(tt[i],a[[j]],b[[j]]))
         }
     }
     lines(tt,f/length(k),col="black",lw=3)
}


get_list_component = function(my_vector, my_index){my_vector[my_index]}
library(RColorBrewer)
mellow.color.pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))

plot_parameters=function(post_sample, par_name){
     possible_num_components = max(lengths(post_sample))
	my_col = mellow.color.pal(possible_num_components)
    
    plot(1:length(post_sample),1:length(post_sample), type="n",main=paste("Posterior Samples of", par_name), ylab=paste(par_name) ,xlab =  "Iterations", ylim = range(unlist(post_sample)))

    
     for(i in 1:possible_num_components){
        which_list_i = which(lengths(post_sample) >= i)
        y = unlist(lapply(post_sample[which_list_i], get_list_component, my_index = i))
        lines(1:length(y), y, col = my_col[i], lwd = 2)
     }
     legend("top", paste("k = ", 1:possible_num_components), col = my_col, lwd = 2, ncol = ceiling(possible_num_components/2))
     
}



