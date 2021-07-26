function [XX, XY] = get_XX_XY(indep_vars, dep_var)
%GET_XX Returns the symbolic (observations as intervals) 
%matrices X'X and X'Y required for the OLS estimator
%   Needs as input a matrix holding all indpendent variables with their 
%   lower and upper bound, respectively, and with only ones in the first 
%   two columns AND a matrix holding the corresponding upper and lower 
%   bounds of the dependent variable for the same set of valid observations.
%   Note: unlike in Billard and Diday (2000), X'X and X'Y are corrected
%   by dividing by m in order to meet the symbolic cross-product
%   expressions derived there. This obviously doesn't change betahat.

m = size(indep_vars,1); % number of valid observations

% a set which contains the column indx of every lower bound variable:
regs_indxset = 1:2:size(indep_vars,2); 

%% X'X %%
XX = zeros(max(regs_indxset), max(regs_indxset));

for j = regs_indxset
    for k = regs_indxset
        for u = 1:m % sum along the u's. Nested loops are executed first
             XX(j,k) = (indep_vars(u,j+1) + indep_vars(u,j)) *...
                 (indep_vars(u,k+1) + indep_vars(u,k)) + XX(j,k);
        end
    end
end
XX = 0.25 / m * XX; 
% Problem: Every even row indx and every even col indx contains only zeros
XX = XX(any(XX,2),:); % delete the zero rows
XX = XX(:,any(XX,1)); % delete the zero cols
% XX is done 

%% X'Y %%
XY = zeros(max(regs_indxset), 1);

for j = regs_indxset
    for u = 1:m
        XY(j,1) = (indep_vars(u,j+1) + indep_vars(u,j)) *...
            (dep_var(u,2) + dep_var(u,1)) + XY(j,1);
    end
end
XY = 0.25 / m * XY;
% Problem: Every even row indx contains only zeros
XY = XY(any(XY,2),:); % delete the zero rows
% XY is done

end

