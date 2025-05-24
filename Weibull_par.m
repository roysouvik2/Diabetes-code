% Function to compute the samples of each weight using Weibull CDF

function val = Weibull_par(N,scale,shape)


% Generate a random permutation of N integers
permut = randperm(N);


% Dividing the range of the Weibull CDF into evenly partitioned subintervals
% to generate probabilities
ran_num = unifrnd((permut-1)/N, permut/N);


% For each probability in the subinterval, find the corresponding value of the weights
% by inverting the Weibull CDF with scale and shape determined above
val = wblinv(ran_num,scale,shape);











