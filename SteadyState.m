function F = SteadyState(L,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s,ME_sy_cd_z, ME_sm_cd_z)
% SteadyState calculates the steady state

F(1) = Lin_Int((L(1)/(L(1)+L(2))), (L(1)+L(2)), ME_sy_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)-L(1);

F(2) = Lin_Int((L(1)/(L(1)+L(2))), (L(1)+L(2)), ME_sm_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)-L(2);


end