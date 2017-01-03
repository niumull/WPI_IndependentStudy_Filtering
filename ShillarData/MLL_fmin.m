function[negLogLike] = MLL_fmin(startvals)
%% 
% Description:
% construct log likelihood function 
% 
% Input: 
% startvals:initial values
% Y: observation value
% X_0:initial guess of X
% Output:
% negLogLike: negative loglikelihod function

%%
% Input the observation value of Y
Y(:,1) = xlsread('Shiller data.xlsx','Sheet3','F3:F66');
Y(:,2) = xlsread('Shiller data.xlsx','Sheet3','E3:E66');
T = length(Y);
X_0 = [-3.2675,0.3570];


% Calculate X series by formula in CIR and S&P model
 X(1,:) = X_0;  % initialization
for i = 2:T+1
   
        X(i,1) = X(i-1,1)+(2*startvals(3)^2*startvals(4)^2-startvals(5)^2)/(2*exp(X(i-1,1)))-startvals(3)^2+...
            startvals(5)/sqrt(exp(X(i-1,1)))*randn(1);
        if X(i,1)<-4.5
            X(i,1) = -4.5;
        else if X(i,1)>-2
                X(i,1) = -2;
            end
        end
        X(i,2) = startvals(6)*exp(X(i-1,1))+startvals(7)*sqrt(exp(X(i-1,1)))*...
            (startvals(8)*randn(1)+sqrt(1-startvals(8)^2)*randn(1));
    
end

% Calculate X_hat = E(Xn|Fn-1);
for j = 1:T+1
        X_hat(j,1) = X(j,1)+(2*startvals(3)^2*startvals(4)^2-startvals(5)^2)/(2*exp(X(j,1)))-startvals(3)^2;
        X_hat(j,2) = startvals(6)*exp(X(j,1));
end

% Construct the loglikelihood function
%X_hat = X_hat(2:end,:);
%X = X(2:end,:);
negLogLike = 0;
for i = 2:T+1
      A = startvals(5)^2/exp(X(i-1,1))+startvals(1)^4;
      B = startvals(7)*startvals(5)*startvals(8);
      C = startvals(7)*startvals(5)*startvals(8);
      D = startvals(7)^2*exp(X(i-1,1))+startvals(2)^4;
      E = A*D-B*C;
      Sigma = [A,B;C,D];
      Sigma_inv = [D/E,-B/E;-C/E,A/E];
      Neu = Y(i-1,:)-X_hat(i,:);
      q = log(det(Sigma))+Neu*Sigma_inv*Neu';
      negLogLike = negLogLike+q;
end
%negLogLike = 0.5*negLogLike;% use 0.5 instead of -0.5 because we want to use fminsearch function
end
