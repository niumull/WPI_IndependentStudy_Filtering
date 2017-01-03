function [negLogLike,d_negLogLike] = MLL_newton(Q1_sqrt,Q2_sqrt,k,theta,sigma,mu,a,rho,L,Y)
%% Description;
% Input: parameters
% Output: negative max likelihood function:negLogLike & differentiate it by
%         each parameter:d_negLogLike

%%
% syms sqrt_Q1 sqrt_Q2 k theta sigma mu a rho 


for i = 1:L
    V(i) = Y(i,1)-Q1_sqrt^2*randn(1);
    A(i) = sigma^2/exp(V(i))+Q1_sqrt^4;
    B(i) = a*sigma*rho;
    C(i) = a*sigma*rho;
    D(i) = a^2*exp(V(i))+Q2_sqrt^4;
    E(i) = A(i)*D(i)-B(i)*C(i);
    negLogLike = -0.5*sum(log(a^2*sigma^2+sigma^2/exp(V(i))*Q2_sqrt^4+...
                 Q1_sqrt^4*a^2*exp(V(i))+Q1_sqrt^4*Q2_sqrt^4-a^2*sigma^2*rho^2)...
                  +[Y(i+1,1)-V(i)+(2*k^2*theta^2-sigma^2)/(2*exp(V(i)))-k^2, Y(i+1,2)-mu*exp(V(i))]...
                  *[D(i)/E(i),-B(i)/E(i);-C(i)/E(i),A(i)/E(i)]*...
               [Y(i+1,1)-V(i)+(2*k^2*theta^2-sigma^2)/(2*exp(V(i)))-k^2;Y(i+1,2)-mu*exp(V(i))]);
% vi=yi-noise
end
d_negLogLike = [diff(negLogLike,Q1_sqrt),diff(negLogLike,Q2_sqrt),diff(negLogLike,k),diff(negLogLike,theta),...
                diff(negLogLike,sigma), diff(negLogLike,mu), diff(negLogLike,a), diff(negLogLike,rho)];
% A   sigma^2/exp(Y(i,1))+Q1^2
% B^2, C^2 a^2*sigma^2*rho^2 
% D a^2*exp(Y(i,1))+Q2^2
% E (sigma^2/exp(Y(i,1))+Q1^2)*(a^2*exp(Y(i,1))+Q2^2)-a^2*sigma^2*rho^2

end

