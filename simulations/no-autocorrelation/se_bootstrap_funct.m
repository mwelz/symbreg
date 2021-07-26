function [SEs] = se_bootstrap_funct(indep_vars,dep_var,R)
%BOOTSTRAP_FUNCT Returns an array of the bootstrap estimates for the OLS
%estimates betahat for given y and x in interval data.
%   indep_vars - contains interval data (lower bound and upper bound for
%                every regressor) and the first two columns contain only ones
%   dep_var    - contains lower and upper bound for the dependent variable
%   R          - number of boostrap runs (about 2,000 should be fine to obtain
%                unbiased standard errors, see Fox (2015)) 
%   Note: The function 'get_XX_XY.m' is required here

m             = size(indep_vars,1); % number of observations
no_regs       = size(indep_vars,2) / 2; % number of regressors incl. constant
beta_arr      = zeros(no_regs, R); % every col contains results of one bootstrap run
indep_vars_bs = indep_vars;
dep_var_bs    = dep_var;

for r = 1:R
    % reshuffle all m obervations (i.e. the rows)
    for i = 1:m
        rand_indx          = round(unifrnd(0.5, m + 0.5 - realmin));
        indep_vars_bs(i,:) = indep_vars(rand_indx,:);
        dep_var_bs(i,:)    = dep_var(rand_indx,:);
    end
    
    % Get X'X and X'Y:
    [XX, XY] = get_XX_XY(indep_vars_bs, dep_var_bs);
    % OLS estimate:
    [betahat] = XX \ XY; % faster than inv()
    % fill up the array:
    beta_arr(:,r) = betahat; 
end
beta_arr = beta_arr'; % the betas shall be contained in the columns
SEs      = std(beta_arr); % std() gives std of every column
SEs      = SEs';

end

