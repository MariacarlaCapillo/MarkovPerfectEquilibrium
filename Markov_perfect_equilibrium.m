%% Title

% MATLAB code for computing the policy outcome function and the 
% Markov-perfect equilibrium of a three-period OLG with a pay-as-you-go 
% system.

% Author: Mariacarla Capillo, University of Bern

% Date: 2018-02-01

%% Clear All Windows

clear;                                                      % clear Workspace
clc;                                                        % clear Command Window
close all;                                                  % close Figures


%% Model Parameters

% Base information

version=1;                                                 % Version number
BASE_PATH = '[INSERT PATH TO MATLAB FOLDER]';       % Base path for saving all figures

% Set values for the parameters and the exogenous variables

t = 100;                                                    % Provisional number of periods
A = 1;                                                      % Economy-wide productivity level
alpha = 0.3;                                                % Output elasticity of capital
beta = 0.85;                                                % Psychological discount factor
a_y = 1;                                                    % Labor productivity of young workers
a_m = 1.3;                                                  % Labor productivity of middle-aged workers
omega_m = 1;                                                % Political weight of middle-aged workers with respect to young workers
omega_r = 1.25;                                             % Political weight of retirees with respect to young workers
borr_cons = -inf;                                           % Borrowing constraint of young workers and middle-aged workers

% Sequence of cohort growth rates
%n = 1.26*ones(1,t);
%n = 1.22*ones(1,t);
%n = 1.18*ones(1,t);
%n = 1.14*ones(1,t);
n = 1.10*ones(1,t);


% Construct the sequence of cohort growth size from the sequence of
% cohort growth rates
n_size = zeros(1,t);                                        % Prepare vector
n_size(t) = 0.2;                                            % Set starting value of cohort size
for c=1:t-1
    n_size(t-c)= n_size(t-c+1)*n(t-c);                      % Construct series of cohort growth size
end

% Note that the first entry in the sequences n and n_size stands for the 
% value in the last period. Going backward in time means going forward in
% in the sequence.

% Specify the conditions for convergence of the political policy function

tolerance = 0.001;                                          % Tolerance for convergence of political policy function
starttol = tolerance+1;                                     % Initial value which initiates the loop

% Define the grid for the tax rate and the state variables

tau_min = 0;                                                % Lowest labor income tax rate to consider
tau_max = 1;                                                % Highest labor income tax rate to consider
tau_points = 90;                                            % Grid size
tau_incr = (tau_max-tau_min)/(tau_points-1);                % Step size
tau = tau_min:tau_incr:tau_max;                             % Sequence of labor income taxes to iterate through

dis_s_min = 0;                                              % Lowest value of the percentage, the young own of the sum of the assets (short: distribution), to consider
dis_s_max = 1;                                              % Highest value of distribution to consider
dis_s_points = 40;                                          % Grid size
dis_s_incr = (dis_s_max-dis_s_min)/(dis_s_points-1);        % Step size
dis_s = dis_s_min:dis_s_incr:dis_s_max;                     % Sequence of distribution to consider

sum_s_min = 0.1;                                            % Lowest value of sum of assets (short: sum) to consider
sum_s_max = 3;                                              % Highest value of sum to consider
sum_s_points = 40;                                          % Grid size
sum_s_incr = (sum_s_max-sum_s_min)/(sum_s_points-1);        % Step size
sum_s = sum_s_min:sum_s_incr:sum_s_max;                     % Sequence of sum to consider


% Preallocate the space for the solution matrices

Sol_cd_lp = zeros(dis_s_points, sum_s_points, tau_points, 6);   % Solution matrix with equilibrium values of savings, taxes and prices for the last period
Sol_cd = zeros(dis_s_points, sum_s_points, tau_points, 9,t);    % Solution matrix with equilibrium values of savings, taxes and prices for all previous periods
Welfare_eq_cd = zeros(dis_s_points, sum_s_points, tau_points,t);% Matrix with welfare values for every equilibrium
MaxPolicy_cd = zeros(dis_s_points, sum_s_points, t);            % Matrix with tax rates that yield highest welfare for every combination of states
Markov_cd_lp = zeros(dis_s_points, sum_s_points, 6);            % Solution matrix with Markov equilibrium values of savings, taxes and prices for the last period
Markov_cd = zeros(dis_s_points, sum_s_points, 9, t);            % Solution matrix with Markov equilibrium values of savings, taxes and prices for all previous periods
Wealth_dis_cd = zeros(dis_s_points,sum_s_points, t);            % Matrix with wealth distribution for every combination of states
Diff_pol_fun_cd = starttol*ones(1,t);                           % Matrix with difference in policy function (needed to compute convergence)
Tau_T=zeros(dis_s_points,sum_s_points);                         % Matrix for convergence purposes
calibration = zeros(dis_s_points,sum_s_points);                 % Matrix for calibration purposes
Reg_coeff=zeros(5,t);

% Set the options for the solver 'fsolve'
options = optimoptions('fsolve', 'Display', 'final');
rng default

% Educated guesses for the inital values of the tax rates which we need for
% the solver 'fsolve'
tau_guess = 0.1;                                            % Guess for the tax rate tomorrow
taut_guess = 0.1;                                           % Guess for the tax rate the day after tomorrow
r_guess = 1.5;                                              % Guess for the interest rate tomorrow
rt_guess = 1.5;                                             % Guess for the interest rate the day after tomorrow
w_guess = 0.3;                                              % Guess for the wage tomorrow
wt_guess = 0.3;                                             % Guess for the wage the day after tomorrow
 
% Lower and upper bound for random guesses
tau_guess_min = 0;                                          % Guess for the tax rate tomorrow
tau_guess_max = 0.7;

% Start in the last period of the economy (conceptually this corresponds to
% period T-1 since this is the first period where we have an economic
% equilibrium)
z=1;

%% Compute the policy outcome function

% In the next sections, we calculate the policy outcome fuction
% numerically. This is the detailed procedure:
% First step: We calculate the equilibrium for each combination of 
% taxes and states and compute the welfare for each equilibrium.
% Second step: We extract the tax rate that yields the highest welfare for 
% each combination of states. We then save and plot this tax rate as a 
% function of the states. We smooth the tax rate by running a quadratic
% regression.
% Third step: We use the smoothed tax rate in order to compute the Markov 
% equilibrium. We plot the different variables in Markov equilibrium.
% Fourthly: We check if convergence of the policy outcome function is 
% reached.

%% First step: Compute the political policy function

% Run the loop as long as the covergence criterion is not satisfied
while Diff_pol_fun_cd(z) > tolerance
% Set the counter for the iteration to 1
status = 1;
% Loop through the different tax values and the different state variables
for v=1:sum_s_points
    for x=1:dis_s_points
        for y=1:tau_points
            % Distinguish between the last period and all previous periods
            % because the equilibrium equations are different in the last
            % period.
            if z == 1
                % Set the initial point of the solver
                L0 = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1), tau(y),dis_s(x),sum_s(v),tau_guess,r_guess,w_guess);
                % For every combination of state variables and tax rates, find the equilibrium
                % of savings, taxes and prices and save it in the solution matrix Sol_cd_lp().
                f = @(L) EquilibriumCD_LP(L,borr_cons,A, alpha, a_y,a_m,beta,omega_m,omega_r,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),n_size(z), tau(y),dis_s(x),sum_s(v));
                [Sol_cd_lp(x,v,y,:),fval,exitflag,output] = fsolve(f,L0,options);
                % In case the solver did not work, try again with a new random starting value
                it = 1;
                while CheckResult(output.message)== 1 && it < 30
                    L0 = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),tau(y),dis_s(x),sum_s(v),(tau_guess_min+(tau_guess_max-tau_guess_min).*rand(1,1)),r_guess,w_guess);
                    [Sol_cd_lp(x,v,y,:),fval,exitflag,output] = fsolve(f,L0,options);
                    it = it+1;
                end
                % If the solver still didn't work, write NaN (not a number)
                if CheckResult(output.message)== 1
                    Sol_cd_lp(x,v,y,:)= NaN;
                end
                % Compute Welfare for each equilibrium
                Welfare_eq_cd(x,v,y,z) = Welfare_cd( Sol_cd_lp(x,v,y,1), Sol_cd_lp(x,v,y,2), Sol_cd_lp(x,v,y,3), Sol_cd_lp(x,v,y,4), Sol_cd_lp(x,v,y,5),A, alpha, a_y,a_m,beta,omega_m,omega_r,n(z+2),n(z+1),n(z), tau(y),dis_s(x),sum_s(v));
                % Display the current iteration
                fprintf('Iteration %i/%i\n', status, dis_s_points*sum_s_points*tau_points)
                % Adjust the counter
                status = status + 1;
            else
                % Set the initial point of the solver
                L0 = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),tau(y),dis_s(x),sum_s(v),tau_guess,r_guess,w_guess,n(z-1),taut_guess,wt_guess,rt_guess);
                % For every combination of state variables and taxes, find the equilibrium
                % of savings, taxes and prices and save it in the solution matrix Sol_cd().
                f = @(L) EquilibriumCD_new(L,borr_cons,A,alpha,a_y,a_m,beta,n(z+2),n(z+1),n(z),n(z-1),n_size(z+2),n_size(z+1),n_size(z),tau(y),dis_s(x),sum_s(v),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s,Reg_coeff(:,z-1), ME_tt_cd_z, ME_wt_cd_z, ME_rt_cd_z);
                [Sol_cd(x,v,y,:,z),fval,exitflag,output] = fsolve(f,L0,options);
                Welfare_eq_cd(x,v,y,z) = Welfare_cd( Sol_cd(x,v,y,1,z), Sol_cd(x,v,y,2,z), Sol_cd(x,v,y,3,z), Sol_cd(x,v,y,5,z), Sol_cd(x,v,y,7,z), A, alpha, a_y,a_m,beta,omega_m,omega_r,n(z+2),n(z+1),n(z), tau(y),dis_s(x),sum_s(v),Sol_cd(x,v,y,4,z), Sol_cd(x,v,y,6,z), Sol_cd(x,v,y,8,z),n(z-1),ME_sm_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s );
                 % In case the solver did not work or the solution doesn't make sense,
                 % try again with a new random starting value
                it = 1;
                while (CheckResult(output.message)== 1 || isinf(Welfare_eq_cd(x,v,y,z)) == 1 || Sol_cd(x,v,y,1,z)>10 || Sol_cd(x,v,y,3,z)>1 || Sol_cd(x,v,y,4,z)>1) && it < 30 
                    L0 = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),tau(y),dis_s(x),sum_s(v),(tau_guess_min+(tau_guess_max-tau_guess_min).*rand(1,1)),r_guess,w_guess,n(z-1),taut_guess,wt_guess,rt_guess);
                    [Sol_cd(x,v,y,:,z),fval,exitflag,output] = fsolve(f,L0,options);
                    Welfare_eq_cd(x,v,y,z) = Welfare_cd( Sol_cd(x,v,y,1,z), Sol_cd(x,v,y,2,z), Sol_cd(x,v,y,3,z), Sol_cd(x,v,y,5,z), Sol_cd(x,v,y,7,z), A, alpha, a_y,a_m,beta,omega_m,omega_r,n(z+2),n(z+1),n(z), tau(y),dis_s(x),sum_s(v),Sol_cd(x,v,y,4,z), Sol_cd(x,v,y,6,z), Sol_cd(x,v,y,8,z),n(z-1),ME_sm_cd_z,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s );
                    it = it+1;
                end
                % If the solver still didn't work or the solution doesn't make sense,
                % write NaN (not a number)
                if (CheckResult(output.message)== 1 || isinf(Welfare_eq_cd(x,v,y,z)) == 1 || Sol_cd(x,v,y,1,z)>10 || Sol_cd(x,v,y,3,z)>1 || Sol_cd(x,v,y,4,z)>1)
                    Sol_cd(x,v,y,:,z) = NaN;
                    Welfare_eq_cd(x,v,y,z) = -Inf;
                end
                % Display the current iteration and its computing time
                fprintf('Iteration %i/%i\n', status, dis_s_points*sum_s_points*tau_points)
                % Adjust the counter
                status = status + 1;
            end
        end
    end
end


%% Second step: Extract Optimal Tax Rate
% Save the tax rates that yield the highest welfare, for each combination
% of state variables

% Extract the welfare from the current period
Welfare_z_cd = Welfare_eq_cd(:,:,:,z);

% Extract the highest welfare and the position of tax rate that yielded it
[MaxWelf_cd, Position_cd]= max(Welfare_z_cd, [], 3, 'omitnan');

% Transform position into actual tax rate
MaxPolicy_cd(:,:,z) = tau(Position_cd);

% Extract this periods policy function in order to plot it
MaxPolicy_z_cd = MaxPolicy_cd(:,:,z);

% Plot the optimal political policy variable as a function of the state
% variables (last periods taxes and the savings of last period's
% middle-aged households)
MaxPol_fig_cd = figure( 'Name', 'The welfare maximizing tax rate' );
surf(sum_s, dis_s, MaxPolicy_z_cd)
% Label axes
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Political Policy Function'
grid on
colorbar
% Save plot in every loop
path = strcat(BASE_PATH, 'Max_Pol_Fig\');
saveas(MaxPol_fig_cd,[path,sprintf('MaxPol_fig_cd_Period_z_%d_%d.png',z,version)]);

%% Smoth Political Policy Function

% Construct the variables in order to run the regression
% Tax_opt = coeff_1*dis_s + coeff_2*dis_s^2+ coeff_3*sum_s + coeff_4*sum_s^2+
% coeff_5*dis_s*sum_s.

% Prepare the dependent variable
MaxPolicy_reg = MaxPolicy_z_cd(:);

% Prepare the independent variables

% Prepare 'dis_s'
dis_reg = repmat(dis_s',dis_s_points);
dis_reg_def = dis_reg(:,1);
% Prepare 'dis_s^2'
dis_reg_def2 = dis_reg_def.^2;
% Prepare 'sum_s'
sum_reg = repmat(sum_s,sum_s_points);
sum_reg_def_2 = sum_reg(:);
sum_reg_def = sum_reg_def_2(1:dis_s_points*sum_s_points);
% Prepare 'sum_s^2'
sum_reg_def2 = sum_reg_def.^2;
% Prepare the interaction 'dis_s*sum_s'
dis_sum_reg = dis_reg_def.*sum_reg_def;
% Run the regression and save the coefficients
x_reg=[dis_reg_def dis_reg_def2 sum_reg_def sum_reg_def2 dis_sum_reg];
Reg_coeff(:,z) = regress(MaxPolicy_reg,x_reg);


%% Third step: Calculate the Markov Equilibrium

% Calculate again the equilibrium for every combination of states, but 
% this time insert the smoothed political policy function we calculated above.

% Set the counter again to one
status = 1;

for v = 1:sum_s_points
    for x = 1:dis_s_points
        % Distinguish between the last period and all previous periods
        % because the equilibrium equations are different in the last
        % period.
        if z == 1
            % Set the initial point of the solver
            L0_m = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v),tau_guess,r_guess,w_guess);
            % Find the equilibrium of savings,prices and taxes and save it in the
            % solution matrix Markov_cd_lp()
            f = @(L) EquilibriumCD_LP(L,borr_cons,A, alpha, a_y,a_m,beta,omega_m,omega_r,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),n_size(z), max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v));
            [Markov_cd_lp(x,v,:),fval,exitflag,output] = fsolve(f,L0_m,options);
            % In case the solver did not work, try again with a new random starting value
            it = 1;
            while CheckResult(output.message)== 1 && it < 30
                L0_m = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v),(tau_guess_min+(tau_guess_max-tau_guess_min).*rand(1,1)),r_guess,w_guess);
                [Markov_cd_lp(x,v,:),fval,exitflag,output] = fsolve(f,L0_m,options);
                it = it+1;
            end
            % If the solver still didn't work, write NaN (not a number)
            if CheckResult(output.message)== 1
                Markov_cd_lp(x,v,:) = NaN;
            end
            % Display the current iteration
            fprintf('ME Iteration %i/%i\n', status, dis_s_points*sum_s_points)
            % Adjust the counter
            status = status + 1;
            % Fill missing values (NaN's). If we do not do this, the linear
            % interpolation in the next step doesn't work.
            Markov_cd_lp = fillmissing(Markov_cd_lp, 'linear');
        else
            % Set the initial point of the solver
            L0_m = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v),tau_guess,r_guess,w_guess,n(z-1),taut_guess,wt_guess,rt_guess);
            % Find the equilibrium of savings,prices and taxes and save it in the
            % solution matrix Markov_cd()
            f = @(L) EquilibriumCD_new(L,borr_cons,A,alpha,a_y,a_m,beta,n(z+2),n(z+1),n(z),n(z-1),n_size(z+2),n_size(z+1),n_size(z),max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s,Reg_coeff(:,z-1), ME_tt_cd_z, ME_wt_cd_z, ME_rt_cd_z);
            [Markov_cd(x,v,:,z),fval,exitflag,output] = fsolve(f,L0_m,options);
            % In case the solver did not work or the solution doesn't make sense,
            % try again with a new random starting value
            it = 1;
            while (CheckResult(output.message)== 1 || Markov_cd(x,v,1,z)>10 || Markov_cd(x,v,3,z)>1 || Markov_cd(x,v,4,z)>1) && it < 30
                L0_m = startingvalue_cd( borr_cons,A, alpha, a_y, a_m, beta,n(z+2),n(z+1),n(z),n_size(z+2),n_size(z+1),max(0,Reg_coeff(1,z)*dis_s(x)+Reg_coeff(2,z)*(dis_s(x))^2+Reg_coeff(3,z)*sum_s(v)+Reg_coeff(4,z)*(sum_s(v))^2+Reg_coeff(5,z)*dis_s(x)*sum_s(v)),dis_s(x),sum_s(v),(tau_guess_min+(tau_guess_max-tau_guess_min).*rand(1,1)),r_guess,w_guess,n(z-1),taut_guess,wt_guess,rt_guess);
                [Markov_cd(x,v,:,z),fval,exitflag,output] = fsolve(f,L0_m,options);
                it = it+1;
            end
            % If the solver still didn't work or the solution doesn't make sense,
            % write NaN (not a number)
            if CheckResult(output.message)== 1 || Markov_cd(x,v,3,z)>1 || Markov_cd(x,v,4,z)>1 || Markov_cd(x,v,1,z)>10
                Markov_cd(x,v,:,z) = NaN;
            end
            % Display the current iteration and its computation time
            fprintf('ME Iteration %i/%i\n', status, dis_s_points*sum_s_points)
            % Adjust the counter
            status = status + 1;
            % Fill missing values (NaN's). If we do not do this, the linear
            % interpolation in the next step doesn't work.
            Markov_cd(:,:,:,z) = fillmissing(Markov_cd(:,:,:,z), 'linear');
        end
    end
end


%% Plot the savings and tax functions in Markov equilibrium

% Plot the savings of the young in Markov equilibrium
if z == 1
    ME_sy_cd_z = Markov_cd_lp(:,:,1);
else
    ME_sy_cd_z = Markov_cd(:,:,1,z);
end
Sav_y_fig_cd = figure( 'Name', 'The savings of the young in Markov equilibrium' );
surf(sum_s, dis_s, ME_sy_cd_z)
% Label axes
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Savings of the young in M.E.'
grid on
colorbar
% Save plot in every loop
path =  strcat(BASE_PATH, 'Saving_young_ME_Fig\');
saveas(Sav_y_fig_cd,[path,sprintf('Sav_y_cd_Period_z_%d_%d.png',z,version)]);

% Plot the savings of the middle-aged agents in Markov equilibrium
if z == 1
    ME_sm_cd_z = Markov_cd_lp(:,:,2);
else
    ME_sm_cd_z = Markov_cd(:,:,2,z);
end
Sav_ma_fig_cd = figure( 'Name', 'The savings of the middle-aged agents in Markov equilibrium' );
surf(sum_s, dis_s, ME_sm_cd_z)
% Label axes
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Savings of the middle-aged workers in M.E.'
grid on
colorbar
% Save plot in every loop
path =  strcat(BASE_PATH, 'Savings_ma_ME_Fig\');
saveas(Sav_ma_fig_cd,[path,sprintf('Sav_ma_cd_Period_z_%d_%d.png',z,version)]);

% Save the tax rate tomorrow in Markov equilibrium
if z ==1
    ME_tt_cd_z = Markov_cd_lp(:,:,3);
else
    ME_tt_cd_z = Markov_cd(:,:,3,z);
end
figure( 'Name', 'The tax tomorrow in Markov equilibrium' );
surf(sum_s, dis_s, ME_tt_cd_z)

% Save the equilibrium wage tomorrow
if z ==1
    ME_wt_cd_z = Markov_cd_lp(:,:,4);
else
    ME_wt_cd_z = Markov_cd(:,:,5,z);
end
figure( 'Name', 'The wage tomorrow in Markov equilibrium' );
surf(sum_s, dis_s, ME_wt_cd_z)
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Wage'
grid on
colorbar

% Save the equilibrium interest tomorrow
if z ==1
    ME_rt_cd_z = Markov_cd_lp(:,:,5);
else
    ME_rt_cd_z = Markov_cd(:,:,7,z);
end
figure( 'Name', 'The interest tomorrow in Markov equilibrium' );
surf(sum_s, dis_s, ME_rt_cd_z)
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Interest rate'
grid on
colorbar

%% Fourthly: Check if convergence is reached
% Check if the regression coefficients stay constant over time If they
% stay constant during 3 periods, we say that the policy outcome function
% converged.

if z<3
    Diff_pol_fun_cd(z+1) = starttol;
else
    %Diff_pol_fun_cd(z+1) = sum(sum(abs(MaxPolicy_cd(:,:,z)-MaxPolicy_cd(:,:,z-1))))+sum(sum(abs(MaxPolicy_cd(:,:,z-1)-MaxPolicy_cd(:,:,z-2))));
    Diff_pol_fun_cd(z+1) = sum(abs(Reg_coeff(:,z)-Reg_coeff(:,z-1)))+sum(abs(Reg_coeff(:,z-1)-Reg_coeff(:,z-2)));
end

% Go back one period
z = z+1;

end

%% Compute steady state

% Use the function 'SteadyState' to find the steady state
fun = @(L) SteadyState(L,dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s,ME_sy_cd_z, ME_sm_cd_z);
% Set the initial value of the solver
x0=[0.2,0.6];
% Solve for the steady state
result = fsolve(fun,x0);

% Display the steady state tax rate
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,3,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,4,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)
% Display the steady state wage
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,5,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,6,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)
% Display the steady state interest rate
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,7,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)
Lin_Int((result(1)/(result(1)+result(2))), (result(1)+result(2)), Markov_cd(:,:,8,z-1),dis_s_points,dis_s_incr,dis_s,sum_s_points,sum_s_incr,sum_s)

%% Plot the smoothed policy outcome function

% Construct the policy outcome for every combination of states by
% multiplying every combination of states with the regression coefficients
SmoothTax=zeros(dis_s_points,sum_s_points);

for x=1:dis_s_points
    for v=1:sum_s_points
        SmoothTax(x,v)=max(0,Reg_coeff(1,z-1)*dis_s(x)+Reg_coeff(2,z-1)*(dis_s(x))^2+Reg_coeff(3,z-1)*sum_s(v)+Reg_coeff(4,z-1)*(sum_s(v))^2+Reg_coeff(5,z-1)*dis_s(x)*sum_s(v));
    end
end

% Plot the policy outcome function
figure( 'Name', 'The policy outcome function' );
surf(sum_s, dis_s, SmoothTax)
% Label axes
xlabel 'Sum of savings'
ylabel 'Distribution'
zlabel 'Policy outcome function'
grid on
colorbar       



% The end
