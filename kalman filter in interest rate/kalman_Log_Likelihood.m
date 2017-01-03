%===========================================================================================================================
%***DESCRIPTION
% This function is to calculate the negative value of log likelihood.  
% Kalman filtering-based calibration (Date & Wang, Ch 3)
% Built a function of Maxlikelihood in order to be used in "fminsearch" function 
% The main function 
%   Y = a_0 + a_1*x_1 + a_2*x_2                                                                       % Date & Wang (3)
%   x_1(n) = -alpha_1/(1+alpha_1)*x_1(n-1) + [sigma_1_1/1+alpha_1, sigma_1_2/1+alpha_1]*[B_1, B_2]'   % Date & Wang (2)
%          = A_1*x_1(n-1)+Q_1*[B_1, B_2]'
%   x_2(n) = -alpha_2/(1+alpha_2)*x_2(n-1) + [sigma_2_1/1+alpha_2, sigma_2_2/1+alpha_2]*[B_1, B_2]'   % Date & Wang (2)
%          = A_2*x_2(n-1)+Q_2*[B_1, B_2]'
%
%***INPUT
%  v:   10*1 matrix contain sigma, alpha, lambda, R and mu
%  y:   data we use
%  x0:  initial guess of estimator
%  P0:  initial guess of variance of x
%  tau: maturity date of a specific bond
% 
%***OUTPUT
% The negative value of log likelihood function
%
%***REFERENCE
% Forecasting, structural time series models and the Kalman filter ANDREW C. HARVEY, 2001 
% Linear Gaussian affine term structure models with unobservable factors: Calibration and yield forecasting, Date&Wang(2008)
% http://faculty.washington.edu/eeholmes/Files/Intro_to_kalman.pdf
%% 
%----------------------------------------------------------------------------------------------------------------------------
function negLogLike = kalman_Log_Likelihood(v0,y,x0,P0,tau)

%Initialize parameters in kalman filter 
sigma_1_1 = sqrt(exp(v0(1)));
sigma_1_2 = sqrt(exp(v0(2)));
sigma_2_1 = sqrt(exp(v0(3)));
sigma_2_2 = sqrt(exp(v0(4)));

alpha_1 = v0(5); 
alpha_2 = v0(6); 

lambda_1 = v0(7);
lambda_2 = v0(8);

R = sqrt(exp(v0(9)));

mu = v0(10);

% Date & Wang equation (2) in discrete time using Backward Euler method 
% Represented by equation X_n+1 = A*X_n + Q*B (B is a standard BM)
Q_1 = [sigma_1_1/(1+alpha_1), sigma_1_2/(1+alpha_1)];
Q_2 = [sigma_2_1/(1+alpha_2), sigma_2_2/(1+alpha_2)];

A_1 = (-alpha_1)/(1+alpha_1);
A_2 = (-alpha_2)/(1+alpha_2);

% Two-factor model's coefficients and constant term(Date & Wang equation (3))a_0, a_1, a_2 in Y
H = @(x) (1-exp(-x))/x;
r_inf = mu + lambda_1*(sigma_1_1/alpha_1+sigma_2_1/alpha_2) + lambda_2*(sigma_1_2/alpha_1+sigma_2_2/alpha_2); % Date & Wang (4)
omega = H(alpha_1*tau)*(lambda_1*sigma_1_1/alpha_1+lambda_2*sigma_1_2/alpha_2-(sigma_1_1^2+sigma_1_2^2)/(alpha_1^2)...
                        -(sigma_2_1*sigma_1_1+sigma_1_2*sigma_2_2)/(alpha_1*alpha_2))...
       +H(alpha_2*tau)*(lambda_1*sigma_2_1/alpha_1+lambda_2*sigma_2_2/alpha_2-(sigma_2_1^2+sigma_2_2^2)/(alpha_2^2)...
                        -(sigma_2_1*sigma_1_1+sigma_1_2*sigma_2_2)/(alpha_1*alpha_2))...
  +0.5*(H(2*alpha_1*tau)*(sigma_1_1^2+sigma_1_2^2)/(alpha_1^2)...
       +H(2*alpha_2*tau)*(sigma_2_1^2+sigma_2_2^2)/(alpha_2^2)...
       +H((alpha_1+alpha_2)*tau)*((2*sigma_2_1*sigma_1_1+2*sigma_1_2*sigma_2_2)/(alpha_1*alpha_2))); % Date & Wang (5)

a_0 = r_inf - omega;
a_1 = -H(alpha_1*tau);
a_2 = -H(alpha_2*tau);


% Initialize 
n=length(y);
x_1=zeros(1,n); P_1=zeros(1,n); x_1_0=zeros(1,n); P_1_0=zeros(1,n);
x_2=zeros(1,n); P_2=zeros(1,n); x_2_0=zeros(1,n); P_2_0=zeros(1,n);
Ft_1=zeros(1,n); vt=zeros(1,n);Ft_2=zeros(1,n);

%Cauculate Ft_1, Ft_2 and vt by kalman filter
for t=1:n
    if t==1
        x_1_0(1) = x0(1); %denotes x_1^0
        x_2_0(1) = x0(2);
        P_1_0(1) = P0(1); %denotes V_1^0
        P_2_0(1) = P0(2); 
    else
        x_1_0(t) = A_1*x_1(t-1); %Harvey 3.2.2a
        x_2_0(t) = A_2*x_2(t-1);
        P_1_0(t) = A_1^2*P_1(t-1) + Q_1*Q_1'; %Harvey 3.2.2b
        P_2_0(t) = A_2^2*P_2(t-1) + Q_2*Q_2';
    end
    Kt_1 = P_1_0(t)*a_1/(a_1^2*P_1_0(t)+R*R);
    Kt_2 = P_2_0(t)*a_2/(a_2^2*P_2_0(t)+R*R);
    Ft_1(t) = a_1^2*P_1_0(t)+R*R;
    Ft_2(t) = a_2^2*P_2_0(t)+R*R;
    vt(t) = y(t)-a_0-a_1*x_1_0(t)-a_2*x_2_0(t);
    x_1(t) = x_1_0(t) + Kt_1*vt(t); %Harvey 3.2.3a
    x_2(t) = x_2_0(t) + Kt_2*vt(t);
    P_1(t) = P_1_0(t)-Kt_1*P_1_0(t)*a_1; %Harvey 3.3.3b
    P_2(t) = P_2_0(t)-Kt_2*P_2_0(t)*a_2;
end

%Calculate negative log likelihood value
negLogLike = (1/2)*sum(vt.^2./Ft_1) + (1/2)*sum(log(abs(Ft_1))) + (n/2)*log(2*pi)...
             +(1/2)*sum(vt.^2./Ft_2) + (1/2)*sum(log(abs(Ft_2))) ;
end
%----------------------------------------------------------------------------------------------------------------------------