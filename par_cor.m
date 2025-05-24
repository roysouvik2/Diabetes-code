% Computing the PRCC of the output of interest and 
% the sensitive weights

function [] = par_cor(samp,num,G_sleep,H_sleep,G0,H0,t_span)

% For each column of the LHS matrix, corresponding to a set of weights, solve the optimal control problem
% and compute the output of interest (in this case the integral of u)

for i = 1:num
    Gfinal(i,1) = ODE_solver(samp(:,i),G_sleep,H_sleep,G0,H0,t_span);
end

% Computing the PRCC coefficients

% Transposing the LHS matrix such that each row now corresponds to a set of weights and pasting the 
% column vector of the output of interest for each sample to this transposed matrix
samp = samp';

mat = [samp Gfinal];


% Computing the Spearman PRCC and p values of this joint matrix
[rho,pval] = partialcorr(mat,'type', 'Spearman');

s = size(samp);

% Outputs: Check the last column of the PRCC matrix rho and the table of p values.
tum = array2table(pval(1:s(2),s(2)+1), ...
    'VariableNames',{'T'},...
    'RowNames',{'Rg','a_g','k_g','K_g','E_g','beta_H','H_v','G_A','G_W'});
fprintf('\n\n')
disp('p-values:')
fprintf('\n')
disp(tum)

fprintf('Correlation Matrix:\n\n')
disp(rho)
