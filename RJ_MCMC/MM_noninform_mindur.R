lambda_hat = 5/100
min_dur = 10

flag_sim = 0; flag_data = 1
flag_first = 0

kk = 57
#outlier_index = which.max(abs(Xin))

#PICKS UP ON THE OUTLIERS for 57, 91


if(flag_sim == 1){
	sim_data = c(sapply(c(0,.05,0), rnorm, sd = sqrt(0.001) ,n = 34 ))
	plot(sim_data)
	N = length(sim_data)
}else if (flag_data == 1){
	
	if(flag_first == 1){
		load('/Users/mjensen1/Dropbox/MM_Limitiation_switch/NewData_Oct24_2018/data_scrubbing_code/MM_data_Oct2018.Rdata')
		source_url("https://raw.githubusercontent.com/jensenmelanie/Molecular-Motor-Modeling/master/data_analysis_fncs")
	
		particle_data = all_kin1DDB
		particle_type = "kin1 DDB"
	}
	
	mp_data = particle_data[[kk]]
	t1 = mp_data[,'time.s']- mp_data[1,'time.s']
	x1 = mp_data[,'x.microns']- mp_data[1,'x.microns']
	y1 = mp_data[,'y.microns']- mp_data[1,'y.microns']
	

	proj_coord = projection_function(t1,x1, y1) 
	Xn = proj_coord$Long; filled = Xn_fill_fnc(t1, Xn, cut_off = 10)
	Xin = diff(filled$Xn_fill)
	#Xin[outlier_index] = (Xin[outlier_index-1] + Xin[outlier_index + 1])/2 + rnorm(n =1, mean = 0 , sd = abs(Xin[outlier_index])/10)
	par(mfrow = c(1,2))
	plot(filled$Xn_fill, type = 'l')
	plot(Xin,pch = 16)
	title(kk)
	sim_data = diff(filled$Xn_fill)
	N= length(sim_data)

}


##=============================================
##------------- Functions
##=============================================
proposal_r_fnc = function(N, min_dur, prob_switch){
	j = 1; tau_vec = c(0)
	while(tau_vec[j]+ min_dur <= N){
		new_tau = rgeom(1, prob = prob_switch) + tau_vec[j] + min_dur
		if(new_tau >= (N- min_dur)){break}
		tau_vec = c(tau_vec, new_tau)
		j = j+1
	}
	r_vec = rep(0, N); r_vec[tau_vec[-1]] = 1
	return(r_vec)
}


state_list_fnc = function(y, r_vec){
	N = length(y); tau_r = which(r_vec == 1)
	loop_tau = c(0, tau_r, N)
	state_index = list(); y_state = list()
	for(st in 1:(length(tau_r)+1)){
		state_index[[st]] = (loop_tau[st]+1):loop_tau[st+1]
		y_state[[st]] = y[state_index[[st]]]
	}

	return(list(index = state_index, y_st = y_state))
}# End of fnc


Sr_fnc = function(y, r_vec){
	N = length(y); tau_r = which(r_vec == 1)
	state_list = state_list_fnc(y, r_vec)
	return(sum(sapply(state_list$y_st,function(y_st){sum((y_st - mean(y_st))^2)})))	
}# End of fnc

loglike_r = function(y, r_vec){
	Kr = sum(r_vec)+1; N = length(y)
	tau_r = which(r_vec == 1); n_k = diff(c(0, tau_r, N)) 
	return(-1*sum(log(n_k))+Kr*log(pi)-1*((N-Kr)/2)*log(Sr_fnc(sim_data, r_vec))+lgamma((N-Kr)/2))
}# End of fnc


lprior_r_dur = function(r_vec, prob_switch, min_dur){
	tau_r = c(0,which(r_vec==1),N)
	A = sum(sapply( diff(tau_r)-1-2*(min_dur-1), function(x){max(x,0)}))
	Kr = sum(r_vec)+1
	return((Kr-1)*log(prob_switch) + A*log(1- prob_switch))
}# End of fnc


r_vec012_fnc = function(r_vec, min_dur){
	N = length(r_vec)
	r_filled = rep(0, N)
	tau_vec = which(r_vec == 1)
	for(j in 1:length(tau_vec)){
		r_filled[max(0,(tau_vec[j] - (min_dur-1))): min((tau_vec[j]+ (min_dur-1)),N)] = 2
	}
	r_filled[(N-min_dur + 1):N] = 2; r_filled[1:(min_dur + 1)] = 2;
	r_filled[tau_vec] = 1
 	return(r_filled)
}# End of fnc




##=============================================
##------------- Algorithm parameters
##=============================================

save_num = 20000
r_mat = matrix(nrow = save_num, ncol = N)
iterations = 50000
accept_ind = matrix(nrow = iterations, ncol = 3)

save_it = 1


r_vec= proposal_r_fnc(N, min_dur, lambda_hat)

##=============================================
##------------- Start of the iteration
##=============================================

for( it in 1: iterations){

	if(it > (iterations- save_num +1)){
		save_it = save_it +1
	}




##--------------- A whole new r vector
	
	r_prop = proposal_r_fnc(N, min_dur, lambda_hat)
	
	like_ratio  = loglike_r(sim_data,r_prop)- loglike_r(sim_data,r_vec)
	#independent sampler so prior ratio and the transition ratio (q) cancel
	
	alpha_1 = min(1,exp(like_ratio))
	
	
	if( alpha_1>= 1){
		r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,1] = 1
	}else if(runif(1, min = 0, max = 1) <= alpha_1){
		r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,1] = 1
	}else{
		r_mat[save_it,] = r_vec; r_vec = r_vec; accept_ind[it,1] = 0
	}
	
##--------------- Add or subtract a single switch point

	r_vec012 = r_vec012_fnc(r_vec, min_dur)
	to_sample = c(which(r_vec012 == 1), which(r_vec012 == 0))
	s_bd = sample(to_sample,1)
	r_prop = r_vec; r_prop[s_bd]  = 1-r_vec[s_bd]
	
	q_rprop_r = 1/length(to_sample)
	q_r_rprop = 1/length(c(which(r_vec012_fnc(r_prop, min_dur) == 0), which(r_vec012_fnc(r_prop, min_dur) == 1)))
	
	log_q_ratio = log(q_r_rprop)- log(q_rprop_r)
	log_prior_ratio = lprior_r_dur(r_prop, lambda_hat, min_dur) - lprior_r_dur(r_vec, lambda_hat, min_dur)
	like_ratio  = loglike_r(sim_data,r_prop)- loglike_r(sim_data,r_vec)
	
	alpha_2 = min(1, exp(like_ratio + log_prior_ratio - log_q_ratio))
	
	if(alpha_2 >=1){
		r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,2] = 1
	}else if(runif(1,0,1) <= alpha_2){
		r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,2] = 1
	}else{
		r_mat[save_it,] = r_vec; r_vec = r_vec; accept_ind[it,2] = 0
	}


	
##--------------- Move the position of a single swithc point
	if(sum(r_vec >0)){

		s_index = sample( which(r_vec ==1),1)
		
		r_nos = r_vec; r_nos[s_index] = 0
		r_vec012 = r_vec012_fnc(r_nos, min_dur)
		
		sprime_to_sample = which(r_vec012 == 0); 
		s_prime_index = sample(sprime_to_sample,1)
		q_rprop_r = 1/length(sprime_to_sample) 	# q(r_prop)|r)

		r_prop = r_nos; r_prop[s_prime_index] = 1
		
		# q(r|r_prop)
		r_prop_prime = r_prop; r_prop_prime[s_prime_index] = 0
		rprop_vec012 = r_vec012_fnc(r_prop_prime, min_dur)
		q_r_rprop = 1/length(which(rprop_vec012 == 0))
	
		#computing the ratio
		log_q_ratio = log(q_r_rprop) - log(q_rprop_r)
		log_prior_ratio = lprior_r_dur(r_prop, lambda_hat, min_dur) - lprior_r_dur(r_vec, lambda_hat, min_dur)
		like_ratio  = loglike_r(sim_data,r_prop)- loglike_r(sim_data,r_vec)
		
		alpha_3 = min(1, exp(like_ratio + log_prior_ratio + log_q_ratio))
	

		if(alpha_3 >=1){
			r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,3] = 1
		}else if(runif(1,0,1) <= alpha_3){
				r_mat[save_it,] = r_prop; r_vec = r_prop; accept_ind[it,3] = 1
			}else{
				r_mat[save_it,] = r_vec; r_vec = r_vec; accept_ind[it,3] = 0
			} # end of the test cases
	}else{
		
		r_mat[save_it,] = r_vec; r_vec = r_vec; accept_ind[it,3] = 0

	} # end of option 3
	
	
} # End of iterations

##=============================================
##-------------plotting
##=============================================





k_ests = apply(r_mat, 1, sum)
probk_est = summary(as.factor(sort((k_ests))))/length(k_ests)
k_map = sort(unique(k_ests))[which.max(probk_est)]



print(apply(accept_ind,2, sum)/nrow(accept_ind))



hist_prob_tau = apply(r_mat, 2, sum)/nrow(r_mat)


tau_map = sort(order(hist_prob_tau, decreasing = TRUE)[1:k_map])

#tau_map = sort(order(hist_prob_tau, decreasing = TRUE)[c(1,2,4)])


tau_map_loop = c(0, tau_map, N); st_index = list()
for(st in 1:(length(tau_map)+1)){
	st_index[[st]] = (tau_map_loop[st]+1):(tau_map_loop[st+1])
}

seg_means = sapply(1:length(st_index), function(index){mean(sim_data[st_index[[index]]])})


layout(matrix(c(1, 1,2,3), 2,2, byrow = TRUE))
plot(sim_data, type = 'l', main = "More Bayesian Approach", ylab = 'position')
#abline(v = 34*c(1:2), lty = 1, col = 'dodgerblue3')
for(st in 1:length(seg_means)){
	lines(st_index[[st]], rep(seg_means[st],length(st_index[[st]])), col = 'springgreen3', lwd = 2)
}


plot(sort(unique(k_ests)), probk_est, type = 'h', ylim = c(0,1), xlab = 'Number of switches', xlim = c(0, max(k_ests)), ylab = "Probability")
points(sort(unique(k_ests)), probk_est, pch = 16)
title("Posterior of K")

plot(apply(r_mat[which(k_ests == k_map),],2,sum)/length(which(k_ests == k_map)), type = 's', main = paste("P(r_i = 1| k = ", k_map, ")"), ylab = "Probability", xlab = 'r_i')

