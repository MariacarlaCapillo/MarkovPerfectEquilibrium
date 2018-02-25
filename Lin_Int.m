function MaxPolicy_Interp = Lin_Int(sy, sm, Matrix_Interpol, s_y_points,s_y_incr, s_y,s_m_points,s_m_incr,s_m)
%Linear_Int_eq linearly interpolates data points in 'Matrix_Interpol'

% Linear_Int_eq(sy,sm, Matrix_Interpol) linearly interpolates the data in
% the input matrix "Matrix_interpol". It takes the two variables sy
% and sm and calculates the corresponding data point in "Matrix_Interpol"
% with linear interpolation.
% The procedure is the following: The function will take the inputs sy and
% sm and check where they lie compared to the vector s_y and s_m. It will
% locate their position and compute the wheigths in order to take a wheighted
% average of existing data points (=linear interpolation).

q = 1;
g = 1;

if isnan(sy)== 1 || isnan(sm)== 1
    MaxPolicy_Interp = NaN;
else
    % Find position of sy
        while sy > s_y(q) && q < s_y_points
            q = q+1;
        end
        if  q > 1
            q = q-1;
        end
        % Calculate wheight for sy
        wh_sy = (sy-s_y(q))/s_y_incr;
        % Find position of sm
        while sm > s_m(g) && g < s_m_points
            g = g+1;
        end
        if g >1
            g = g-1;
        end
        % Calculate wheight for sm
        wh_sm = (sm-s_m(g))/s_m_incr;
        % Linear Interpolation
        MaxPolicy_Interp = Matrix_Interpol(q,g)+wh_sy*(Matrix_Interpol(q+1,g)-Matrix_Interpol(q,g))+wh_sm*(Matrix_Interpol(q,g+1)-Matrix_Interpol(q,g));
end

end

