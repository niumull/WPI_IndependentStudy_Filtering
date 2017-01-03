function [Y_estimate,X_estimate,rmse] = test_rolling(par,Y)
%% Description:

%Test the result using test dataset from year 2010-2015
Y0 = zeros(1,2);
Y0(1) = xlsread('Shiller data.xlsx','Sheet3','F66');
Y0(2) = xlsread('Shiller data.xlsx','Sheet3','E66');
Y_ob = [Y0;Y];
L = length(Y_ob);

Y_bar(:,1) = Y_ob(:,1)-par(1)*randn(L,1);
Y_bar(:,2) = Y_ob(:,2)-par(2)*randn(L,1);

for i = 1:L
    X_estimate(i,1) = Y_bar(i,1)+(2*par(3)^2*par(4)^2-par(5)^2)/(2*exp(Y_bar(i,1)))-par(3)^2+...
            par(5)/sqrt(exp(Y_bar(i,1)))*randn(1);
        
        if X_estimate(i,1)<-4.5
            X_estimate(i,1) = -4.5;
        else if X_estimate(i,1)>-10
                X_estimate(i,1) = -10;
            end
        end
        
    X_estimate(i,2) = par(6)*exp(Y_bar(i,1))+par(7)*sqrt(exp(Y_bar(i,1)))*...
            (par(8)*randn(1)+sqrt(1-par(8)^2)*randn(1)); 
end
X_estimate = X_estimate(1:end-1,:);

%Calculate rmse using test dataset from 2010-2015 
for m = 1:2
    Y_estimate(:,m) = X_estimate(:,m)+par(m)*randn(L-1,1);
    rmse(:,m) = sqrt(sum((Y(:,m)-Y_estimate(:,m)).^2)/(L-1));
end



end