function [x] = newton(startvals,tol,nmax,h,Y_obs)
%% Description:
% Newton Method

%%
% Define neg loglikelihood function and its derivatives in negLogLike.m
Y = [[-3.26754,0.35703];Y_obs];	
L = length(Y_obs);

% Input  parameters in function
k = startvals(3);theta=startvals(4);sigma = startvals(5);mu=startvals(6);
a=startvals(7);rho=startvals(8);
Q1_sqrt=startvals(1); Q2_sqrt=startvals(2);
% Y(:,1) = Y(:,1)-Q1_sqrt^2*randn(L+1,1);
% Y(:,2) = Y(:,2)-Q2_sqrt^2*randn(L+1,1);

syms Q1_sqrt Q2_sqrt k theta sigma mu a rho
[negLogLike,d_negLogLike] = MLL_newton(Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho,L,Y);

%%
digits(4);
x = startvals - double(vpa(subs(negLogLike,[Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho],startvals))\...
                       vpa(subs(d_negLogLike,[Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho],startvals)));
ex = sum(abs(x-startvals));

while ex >= tol && h <= nmax
    x_new = x - double(vpa(subs(negLogLike,[Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho],x))\...
                       vpa(subs(d_negLogLike,[Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho],x)));
%     disp(vpa(subs(negLogLike,[Q_1,Q_2,k,theta,sigma,mu,a,rho],x)))
%     disp(x)
    ex = sum(abs(x_new-x));
    h = h+1;
    x = x_new;
    if x(8)>1 || x(8)<-1
        x(8)=-0.84;
    end
    disp(x)
    disp(ex)
    disp(h)
end
x(1)=x(1)^2;
x(2)=x(2)^2;
end

