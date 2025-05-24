function samp = LHS(theta,num)


len = length(theta); % Length of the weight vector

% This step is for each weight
for i = 1:len
    
    % Generating the Weibull CDF scale and shape parameters
    
    % Mean value
    m = theta(i);
    
    % Standard deviation
    stdev = 0.2*m;
    
    % Newton's method to determine the scale and shape of the Weibull CDF
    [scale,shape] = Newton(m,stdev);
    
    % Generating samples of a weight by inverting the Weibull CDF
    samp(i,:) = Weibull_par(num,scale,shape); 
end



