function L0 = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,coh_grow_yes,coh_grow_tod,coh_grow_tom, coh_size_yes, coh_size_tod, tax_tod,dis_s_state,sum_s_state,tau_guess,r_guess,w_guess,coh_grow_datom,taut_guess,wt_guess,rt_guess)
% startingvalue_cd calculates the vector of starting values for the economic
% equilibrium that is solved with fsolve.

% Distinguish between the last period and all previous periods
% because the equilibrium equations are different in the last
% period.
if nargin <= 17
    default_coh_grow_datom = 0;
    default_taut_guess = 0;
    default_wt_guess = 0;
    default_rt_guess = 0;
else
    default_coh_grow_datom = coh_grow_datom;
    default_taut_guess = taut_guess;
    default_wt_guess = wt_guess;
    default_rt_guess = rt_guess;
end

% Use the guesses to compute the starting value for optimal savings
if nargin <= 17
    sy_start = max(borr_cons,((beta*r_guess*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_y*(1-tax_tod)-w_guess*a_m*(1-tau_guess))/((1+beta)*r_guess)));
else
    sy_start = max(borr_cons,((beta+beta^2)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_y*(1-tax_tod))/(1+beta+beta^2)-(w_guess*a_m*(1-tau_guess))/((1+beta+beta^2)*r_guess)-((a_m*coh_grow_tom+a_y*default_coh_grow_datom*coh_grow_tom)*default_wt_guess*default_taut_guess)/((1+beta+beta^2)*r_guess*default_rt_guess));
end
sm_start = max(borr_cons,((beta*r_guess*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_m*(1-tax_tod)+beta*r_guess*(A*alpha*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^(alpha-1))*dis_s_state*sum_s_state-tau_guess*w_guess*(a_m*coh_grow_tod+a_y*coh_grow_tom*coh_grow_tod))/((1+beta)*r_guess)));
tax_start = tau_guess;
w_start = w_guess;
r_start = r_guess;
k_start = max(5,coh_size_tod*sy_start+coh_size_yes*sm_start);
taxt_start = default_taut_guess;
wt_start = default_wt_guess;
rt_start = default_rt_guess;

% Construct the vector with the starting values
if nargin <= 17
    L0 = [sy_start,sm_start,tax_start,w_start,r_start,k_start];
else
    L0 = [sy_start, sm_start,tax_start,taxt_start, w_start, wt_start, r_start, rt_start,k_start];
end

end

