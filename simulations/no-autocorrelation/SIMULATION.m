clear, clc
tic
%%                      symbreg_big_simulation.m                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a program which creates symbolic data (here: data in intervals),
% performs symbolic regression on the created data and then estimates the
% standard errors of the OLS estimates of the symbolic regression through 
% bootstrapping. We can pepeat this process R times.
% The same random numbers in each single simulation run for the data  
% creation are used, but the bootstrapping draws remain completely random.
% The functions 'get_XX_XY.m' and 'se_bootstrap_funct.m' are required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Initialize
N_set         = 100; % [20,50,100]; % number of observations
alpha_set     = [0, 5];
beta_set      = [-2, 0, 2];
sigma_eps_set = [0.5, 2];
sigma_z_set   = [0.5, 2];
sigma_w_set   = [0.5, 2];

B             = 2000; % number of bootstrap runs
R             = 1000; % number of simulation runs

TotalIters    = length(N_set)*length(alpha_set)*length(beta_set)*...
    length(sigma_eps_set)*length(sigma_z_set)*length(sigma_w_set); % total numer of iterations within one simulation run
PrvResults    = zeros(TotalIters, 10);

% Big Simulation begins here
for r = 1:R
    
    disp(r)            % to see in which sim run we are
    counter       = 0; % counts the iterations

    % to store the values of the i-th simulation result, i=1,..,TotalIters:
    N_val         = zeros(TotalIters,1);
    alpha_val     = zeros(TotalIters,1);
    beta_val      = zeros(TotalIters,1);
    sigsq_eps_val = zeros(TotalIters,1);
    sigsq_z_val   = zeros(TotalIters,1);
    sigsq_w_val   = zeros(TotalIters,1);
    beta0_val     = zeros(TotalIters,1);
    beta1_val     = zeros(TotalIters,1);
    SE_beta0_val  = zeros(TotalIters,1);
    SE_beta1_val  = zeros(TotalIters,1);

    %% 2) Begin the simulation
    for i = N_set
        for j = alpha_set
            for k = beta_set
                for l = sigma_eps_set
                    for m = sigma_z_set
                        for n = sigma_w_set
                         
                            % Initialize the simulation run
                            N           = i; 
                            alpha       = j;
                            beta        = k;
                            sigmasq_eps = l;
                            sigmasq_z   = m;
                            sigmasq_w   = n;
                        
                            %% 3) Create the data
                            % Note: all of the created arrays are col vectors
                            rng(r); % use same rands within one simulation run
                            x_arr        = normrnd(0, sqrt(1), [N,1]);           % Nx1, contains the x_i's
                            eps_arr      = normrnd(0, sqrt(sigmasq_eps), [N,1]); % Nx1, contains the epsilon_i's
                            z_arr_lobo   = normrnd(0, sqrt(sigmasq_z), [N,1]);   % Nx1, contains the z_i's for the lower bound
                            z_arr_upbo   = normrnd(0, sqrt(sigmasq_z), [N,1]);   % Nx1, contains the z_i's for the upper bound
                            w_arr_lobo   = normrnd(0, sqrt(sigmasq_w), [N,1]);   % Nx1, contains the w_i's for the lower bound
                            w_arr_upbo   = normrnd(0, sqrt(sigmasq_w), [N,1]);   % Nx1, contains the w_i's for the upper bound

                            % Create the y_i's according to the DGP y_i = alpha + beta*x_i + epsilon_i:
                            y_arr = alpha + beta * x_arr + eps_arr;         % Nx1, contains the y_i's
    
                            %% 4) Make intervals around the x_i's and y_i's
                            % Upper and lower bound for the y_i's:
                            y_arr_lobo = y_arr - abs(z_arr_lobo);
                            y_arr_upbo = y_arr + abs(z_arr_upbo);
                            % Upper and lower bound for the x_i's:
                            x_arr_lobo = x_arr - abs(w_arr_lobo);
                            x_arr_upbo = x_arr + abs(w_arr_upbo);

                            %% 5) Symbolic regression
                            % Preparing for 'get_XX_XY.m':
                            indep_vars   = [ones(N,1), ones(N,1), x_arr_lobo, x_arr_upbo];
                            dep_var      = [y_arr_lobo, y_arr_upbo];
                            % Get betahat for original data:
                            [XX, XY]     = get_XX_XY(indep_vars, dep_var);
                            betahat_orig = XX \ XY; % faster than inv(XX)

                            %% 6) Bootstrapping to obtain SE(betahat_orig)
                            rng('shuffle') % make bootstrapping draws completely random
                            [SEs] = se_bootstrap_funct(indep_vars,dep_var,B);

                            %% 7) Store the results
                            counter = counter + 1;
                            N_val(counter,1)         = N;
                            alpha_val(counter,1)     = alpha;
                            beta_val(counter,1)      = beta;
                            sigsq_eps_val(counter,1) = sigmasq_eps;
                            sigsq_z_val(counter,1)   = sigmasq_z;
                            sigsq_w_val(counter,1)   = sigmasq_w;
                            beta0_val(counter,1)     = betahat_orig(1);
                            beta1_val(counter,1)     = betahat_orig(2);
                            SE_beta0_val(counter,1)  = SEs(1);
                            SE_beta1_val(counter,1)  = SEs(2);
                                                
                        end
                    end
                end
            end
        end
    end

    % Make a matrix with the results of one run:
    NewResults = [N_val,alpha_val,beta_val,sigsq_eps_val,sigsq_z_val,...
        sigsq_w_val,beta0_val,beta1_val,SE_beta0_val,SE_beta1_val];
    % cumulative results across simulation runs (first 6 cols contain only
    % parameter values):
    NewResults(:,7) = NewResults(:,7) + PrvResults(:,7);
    PrvResults = NewResults;
end % r in R

% Make a table with the final results:
PrvResults(:,7) = PrvResults(:,7) / R; % take means across number of simulatioms
dlmwrite('simulation_output.txt', PrvResults,'precision','%.16f'); % store results in a newly created txt file with 16 decimal places precision
FinalResultsTable = array2table(PrvResults);
FinalResultsTable.Properties.VariableNames = {'N','alpha','beta','sigsq_e',...
    'sigsq_z','sigsq_w','beta0','beta1','SE_0','SE_1'};

toc