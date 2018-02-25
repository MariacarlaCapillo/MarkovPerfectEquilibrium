function F = EquilibriumCD_new(L,borr_cons,A,alpha,a_y,a_m,beta,coh_grow_yes,coh_grow_tod,coh_grow_tom,coh_grow_datom,coh_size_yes,coh_size_tod,coh_size_tom,tax_tod,dis_s_state,sum_s_state,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s,Reg_coeff, ME_tt_cd_z, ME_wt_cd_z, ME_rt_cd_z)
%EquilibriumCD calculates the equilibrium between savings, taxes and prices
%starting in z=2

% EquilibriumCD calculates the equilibrium between savings of the young,
% savings of the middle-aged agents, the tax rate tomorrow, the tax rate
% the day after tomorrow, the interest rate tomorrow, the interest rate the day
% after tomorrow, the wage tomorrow as well as the wage the day after 
% tomorrow starting in period T-2 (z=2) of a
% finite-horizon closed economy. It takes as
% given the tax rate today and the state variables (savings of the young
% yesterday and savings of the middle-aged agents yesterday).

% The savings of the young
F(1) = max(borr_cons,((beta+beta^2)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_y*(1-tax_tod))/(1+beta+beta^2)-(L(5)*a_m*(1-L(3)))/((1+beta+beta^2)*L(7))-((a_m*coh_grow_tom+a_y*coh_grow_datom*coh_grow_tom)*L(6)*L(4))/((1+beta+beta^2)*L(7)*L(8)))-L(1);
% The savings of the middle-aged agents
F(2) = max(borr_cons,((beta*L(7)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*a_m*(1-tax_tod)+beta*L(7)*(A*alpha*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^(alpha-1))*dis_s_state*sum_s_state-L(5)*L(3)*(a_m*coh_grow_tod+a_y*coh_grow_tom*coh_grow_tod))/((1+beta)*L(7))))-L(2);
% The tax rate tomorrow
F(3) = max(0,Reg_coeff(1)*L(1)/(L(1)+L(2))+Reg_coeff(2)*(L(1)/(L(1)+L(2)))^2+Reg_coeff(3)*(L(1)+L(2))+Reg_coeff(4)*(L(1)+L(2))^2+Reg_coeff(5)*(L(1)+L(2))*L(1)/(L(1)+L(2)))-L(3);
% The tax rate the day after tomorrow
F(4) = max(0,Lin_Int((L(1)/(L(1)+L(2))), (L(1)+L(2)), ME_tt_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s))-L(4);
% The equilibrium wage tomorrow
F(5) = (A*(1-alpha)*(L(9)/(a_y*coh_size_tom+a_m*coh_size_tod))^alpha)-L(5);
% The equilibrium wage the day after tomorrow
F(6) = Lin_Int((L(1)/(L(1)+L(2))), (L(1)+L(2)), ME_wt_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)-L(6);
% The equilibrium interest rate tomorrow
F(7) = (A*alpha*(L(9)/(a_y*coh_size_tom+a_m*coh_size_tod))^(alpha-1))-L(7);
% The equilibrium interest rate the day after tomorrow
F(8) = Lin_Int((L(1)/(L(1)+L(2))), (L(1)+L(2)), ME_rt_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)-L(8);
% The capital stock tomorrow
F(9) = coh_size_tod*L(1)+coh_size_yes*L(2)-L(9);
end