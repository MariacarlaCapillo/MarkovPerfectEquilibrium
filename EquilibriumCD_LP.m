function F = EquilibriumCD_LP(L,borr_cons,A, alpha, a_y,a_m,beta,omega_m,omega_r,coh_grow_yes,coh_grow_tod,coh_grow_tom,coh_size_yes,coh_size_tod,coh_size_tom, tax_tod,dis_s_state,sum_s_state)
%EquilibriumCD_LP calculates the equilibrium of savings, taxes and prices
%in z=1

% EquilibriumCD_LP is part of the solver 'fsolve' and
% calculates the equilibrium between savings of the young,
% savings of the middle-aged agents, the tax rate tomorrow, the interest rate
% tomorrow as well as the wage tomorrow in period T-1 (z=1) of a
% finite-horizon closed economy. It takes as
% given the tax rate today and the state variables (savings of the young
% yesterday and savings of the middle-aged agents yesterday). The
% period T-1 differs from all the previous periods, since in this period young
% households will anticipate death in the next period and therefore do not take
% retirement into account.

% The savings of the young
F(1) = max(borr_cons,((beta*L(5)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_y*(1-tax_tod)-L(4)*a_m*(1-L(3)))/((1+beta)*L(5))))-L(1);
% The savings of the middle-aged agents
F(2) = max(borr_cons,((beta*L(5)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_m*(1-tax_tod)+beta*L(5)*(A*alpha*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^(alpha-1))*dis_s_state*sum_s_state-L(3)*L(4)*(a_m*coh_grow_tod+a_y*coh_grow_tom*coh_grow_tod))/((1+beta)*L(5))))-L(2);
% The tax rate tomorrow
F(3) = max(0,Tau_period_T_cd(L(1),L(2),L(4),L(5),a_y,a_m,omega_m,omega_r,coh_grow_yes,coh_grow_tod,coh_grow_tom))-L(3);
% The equilibrium wage tomorrow
F(4) = (A*(1-alpha)*(L(6)/(a_y*coh_size_tom+a_m*coh_size_tod))^alpha)-L(4);
% The equilibrium interest rate tomorrow
F(5) = (A*alpha*(L(6)/(a_y*coh_size_tom+a_m*coh_size_tod))^(alpha-1))-L(5);
% The capital stock tomorrow
F(6) = (coh_size_tod*L(1)+coh_size_yes*L(2))-L(6);

end
