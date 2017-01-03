function [Y_estimate,X_estimate,rmse] = test_09(par,Y)
%% Description:

%%
%Test the result using test dataset from year 2010-2015
L = length(Y);
X_0 = zeros(1,2);
X_0(1) = xlsread('Shiller data.xlsx','Sheet3','F66');
X_0(2) = xlsread('Shiller data.xlsx','Sheet3','E66');
X_0(1) = X_0(1)-par(1)*randn(1);
X_0(2) = X_0(2)-par(2)*randn(1);

X_estimate(1,:) = X_0;

for i = 2:L+1
        X_estimate(i,1) = X_estimate(i-1,1)+(2*par(3)^2*par(4)^2-par(5)^2)/(2*exp(X_estimate(i-1,1)))-par(3)^2+...
            par(5)/sqrt(exp(X_estimate(i-1,1)))*randn(1);
        
        if X_estimate(i,1)<-4.5 % 
            X_estimate(i,1) = -4.5;
        else if X_estimate(i,1)>-10
                X_estimate(i,1) = -10;
            end
        end
        X_estimate(i,2) = par(6)*exp(X_estimate(i-1,1))+par(7)*sqrt(exp(X_estimate(i-1,1)))*...
            (par(8)*randn(1)+sqrt(1-par(8)^2)*randn(1));    
end

X_estimate = X_estimate(2:end,:);

%Calculate rmse using test dataset from 2010-2015 
for m = 1:2
    Y_estimate(:,m) = X_estimate(:,m)+par(m)*randn(L,1);
    rmse(:,m) = sqrt(sum((Y(:,m)-Y_estimate(:,m)).^2)/L);
end

% figure(10)
% hold on
% plot([2010:1:2015],Y_estimate(:,1));
% plot([2010:1:2015],Y(:,1));

%ylim([-1 1])
% legend('Estimated Value','Observations');
% title('Comparison Between Estimated Dividend Yield and Observations')
% xlabel('Year')
% ylabel('X')
% grid ;
% hold off
% 
% figure(20)
% hold on
% plot([2010:1:2015],Y_estimate(:,2));
% plot([2010:1:2015],Y(:,2));
% %ylim([-1 1])
% legend('Estimated Value','Observations');
% title('Comparison Between Estimated S&P Returns and Observations')
% xlabel('Year')
% ylabel('Y')
% grid ;
% hold off

end
