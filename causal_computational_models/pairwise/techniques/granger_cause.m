function [F,c_v] = granger_cause(y,x,alpha,max_lag)
% [F,c_v] = granger_cause(x,y,alpha,max_lag)
% Granger Causality test
% Does Y Granger Cause X?
%
% User-Specified Inputs:
%   x -- A column vector of data
%   y -- A column vector of data
%   alpha -- the significance level specified by the user
%   max_lag -- the maximum number of lags to be considered
% User-requested Output:
%   F -- The value of the F-statistic
%   c_v -- The critical value from the F-distribution
%
% The lag length selection is chosen using the Bayesian information
% Criterion 
% Note that if F > c_v we reject the null hypothesis that y does not
% Granger Cause x

% Chandler Lutz, UCR 2009
% Questions/Comments: chandler.lutz@email.ucr.edu
% $Revision: 1.0.0 $  $Date: 09/30/2009 $
% $Revision: 1.0.1 $  $Date: 10/20/2009 $
% $Revision: 1.0.2 $  $Date: 03/18/2009 $

% References:
% [1] Granger, C.W.J., 1969. "Investigating causal relations by econometric
%     models and cross-spectral methods". Econometrica 37 (3), 424–438.

% Acknowledgements:
%   I would like to thank Mads Dyrholm for his helpful comments and
%   suggestions

%Make sure x & y are the same length
if (length(x) ~= length(y))
    error('x and y must be the same length');
end

%Make sure x is a column vector
[a,b] = size(x);
if (b>a)
    %x is a row vector -- fix this
    x = x';
end

%Make sure y is a column vector
[a,b] = size(y);
if (b>a)
    %y is a row vector -- fix this
    y = y';
end



%Make sure max_lag is >= 1
if max_lag < 1
    error('max_lag must be greater than or equal to one');
end

%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of x

T = length(x);
xm = x;
ym = y;
len = 21;
count = idivide(T, int32(21));
c = double(count);
T = len;

BIC = zeros(max_lag,1);

%Specify a matrix for the restricted RSS
RSS_R = zeros(max_lag,1);

i = 1;
while i <= max_lag
    xstar_ = []; ystar_ = [];
    for k=1:count
        x = xm((k-1)*21 + 1:k*21, :);
        y = ym((k-1)*21 + 1:k*21, :);
        ystar = x(i+1:T,:);
        xstar = [ones(T-i,1) zeros(T-i,i)];
        %Populate the xstar matrix with the corresponding vectors of lags
        j = 1;
        while j <= i
            xstar(:,j+1) = x(i+1-j:T-j);
            j = j+1;
        end
        xstar_ = cat(1, xstar_, xstar);
        ystar_ = cat(1, ystar_, ystar);
    end
    xstar = xstar_;
    ystar = ystar_;
    
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [~,~,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
%     BIC(i,:) = (T-i)*c*log((r'*r)/((T-i)*c)) + (i+1)*c*log((T-i)*c);
    BIC(i,:) = T*log(r'*r/(T*c)) + (i+1)*log(T*c);
    %Put the restricted residual sum of squares in the RSS_R vector
    RSS_R(i,:) = r'*r;
    
    i = i+1;    
end

[~,x_lag] = min(BIC);

%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of y

BIC = zeros(max_lag,1);

%Specify a matrix for the unrestricted RSS
RSS_U = zeros(max_lag,1);

i = 1;
while i <= max_lag
    xstar_ = []; ystar_ = [];
    for k=1:count
        x = xm((k-1)*21 + 1:k*21, :);
        y = ym((k-1)*21 + 1:k*21, :);
        ystar = x(i+x_lag+1:T,:);
        xstar = [ones(T-(i+x_lag),1) zeros(T-(i+x_lag),x_lag+i)];
        %Populate the xstar matrix with the corresponding vectors of lags of x
        j = 1;
        while j <= x_lag
            xstar(:,j+1) = x(i+x_lag+1-j:T-j,:);
            j = j+1;
        end
        %Populate the xstar matrix with the corresponding vectors of lags of y
        j = 1;
        while j <= i
            xstar(:,x_lag+j+1) = y(i+x_lag+1-j:T-j,:);
            j = j+1;
        end
        xstar_ = cat(1, xstar_, xstar);
        ystar_ = cat(1, ystar_, ystar);
    end
    xstar = xstar_;
    ystar = ystar_;
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [~,~,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
    BIC(i,:) = T*log(r'*r/(T*c)) + (i+1)*log(T*c);
%     BIC(i,:) = (T-i-x_lag)*c*log((r'*r)/((T-i-x_lag)*c)) + (i+x_lag+1)*c*log((T-i-x_lag)*c);
    RSS_U(i,:) = r'*r;
    
    i = i+1;
end

[~,y_lag] = min(BIC);

%The numerator of the F-statistic
F_num = ((RSS_R(x_lag,:) - RSS_U(y_lag,:))/y_lag);

%The denominator of the F-statistic
F_den = RSS_U(y_lag,:)/(T-(x_lag+y_lag+1));

%The F-Statistic
F = F_num/F_den;

c_v = finv(1-alpha,y_lag,(T-(x_lag+y_lag+1)));