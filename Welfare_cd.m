function Welf_cd = Welfare_cd( sy, sm, tax_tom, wage_tom, interest_tom, A, alpha, a_y,a_m,beta,omega_m,omega_r,coh_grow_yes,coh_grow_tod,coh_grow_tom, tax_tod,dis_s_state,sum_s_state,tax_datom, wage_datom, interest_datom,coh_grow_datom,ME_sm_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s )
%Welfare_cd computes the Welfare of a given allocation in a given period

% Distinguish between the last period and all previous periods
% because the welfare function is different in the last period.
if nargin > 18
    defaultTaxDatom = tax_datom;
    defaultWageDatom = wage_datom;
    defaultInterestDatom = interest_datom;
    defaultCohgrowDatom = coh_grow_datom;
    defaultMEsm = ME_sm_cd_z;
    defaultDisspoints = dis_s_points;
    defaultDissincr = dis_s_incr;
    defaultDiss = dis_s;
    defaultSumspoints = sum_s_points;
    defaultSumsincr = sum_s_incr;
    defaultSums = sum_s;
else
    defaultTaxDatom = 0;
    defaultWageDatom = 0;
    defaultInterestDatom = 0;
    defaultCohgrowDatom = 0;
    defaultMEsm = 0;
    defaultDisspoints = 0;
    defaultDissincr = 0;
    defaultDiss = 0;
    defaultSumspoints = 0;
    defaultSumsincr = 0;
    defaultSums = 0;
end
     
% The Utility of a young household
c_y_y = a_y*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*(1-tax_tod)-sy;

if nargin <= 18
    c_m_y = a_m*wage_tom*(1-tax_tom)+interest_tom*sy;
else
    c_m_y = a_m*wage_tom*(1-tax_tom)+interest_tom*sy - Lin_Int((sy/(sy+sm)), (sy+sm), defaultMEsm,defaultDisspoints,defaultDissincr,defaultDiss,defaultSumspoints,defaultSumsincr,defaultSums);
end

if nargin > 18
    c_r_y = (a_m*coh_grow_tom+a_y*defaultCohgrowDatom*coh_grow_tom)*defaultWageDatom*defaultTaxDatom + defaultInterestDatom*Lin_Int((sy/(sy+sm)), (sy+sm), defaultMEsm,defaultDisspoints,defaultDissincr,defaultDiss,defaultSumspoints,defaultSumsincr,defaultSums);
end

% If at least one consumption level is negative, set the utility to -inf.
if nargin <= 18
    if c_y_y <= 0 || c_m_y <= 0
        U_y = -inf;
    else
    U_y = log(c_y_y)+beta*log(c_m_y);
    end
else
    if c_y_y <= 0 || c_m_y <= 0 || c_r_y <= 0
        U_y = -inf;
    else
    U_y = log(c_y_y)+beta*log(c_m_y)+ beta^2*log(c_r_y);
    end
end

% The Utility of a middle-aged household
c_m_m = a_m*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*(1-tax_tod)+(A*alpha*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^(alpha-1))*dis_s_state*sum_s_state-sm;
c_r_m = (a_m*coh_grow_tod+a_y*coh_grow_tom*coh_grow_tod)*wage_tom*tax_tom+interest_tom*sm;

% If at least one consumption level is negative, set the utility to -inf.
if c_m_m <= 0 || c_r_m <= 0
    U_m = -inf;
else
    U_m = log(c_m_m)+beta*log(c_r_m);
end

% The Utility of a retiree
c_r_r = (A*alpha*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^(alpha-1))*(1-dis_s_state)*sum_s_state+(a_m*coh_grow_yes+a_y*coh_grow_tod*coh_grow_yes)*(A*(1-alpha)*((1/(a_y*coh_grow_tod+a_m))*(dis_s_state*sum_s_state+((1-dis_s_state)*sum_s_state)/coh_grow_yes))^alpha)*tax_tod;

% If at least one consumption level is negative, set the utility to -inf.
if c_r_r <= 0
    U_r = -inf;
else
    U_r = log(c_r_r);
end

%The Welfare

Welf_cd = U_y+(omega_m/coh_grow_tod)*U_m +(omega_r/(coh_grow_tod*coh_grow_yes))*U_r;

end

