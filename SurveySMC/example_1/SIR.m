
function [x_n,w,x_T] = SIR(Y_n,x_initial,N,T,mu,phi,sigma2)
x_n = zeros(N,T); logWeight = zeros(N,T); weights = zeros(N,T);
N_eff = zeros(1,T); x_n_new = zeros(N,T);
% LogLike = 0;  Loglikelyhood_1 = zeros(1,T_n);
x_n(:,1) = x_initial;
weights(:,1) = (1/N)*ones(N,1);
for t = 2:T
    x_n(:,t) = mu+phi*(x_n(:,t-1)-mu)+sqrt(sigma2)*randn(N,1);
    logWeight(:,t) = -0.5*log(2*pi)-0.5*x_n(:,t)-0.5*(1./exp(x_n(:,t)))*(Y_n(t)^2)+ log(weights(:,t-1));
    weights(:,t) = exp(logWeight(:,t)-max(logWeight(:,t)));
    weights(:,t) = weights(:,t)./sum(weights(:,t));
    N_eff(t) = 1/(sum(weights(:,t).^2));
    if N_eff < N
        [x_n(:,t),weights(:,t),~] = resample(x_n(:,t),weights(:,t),N);
    end
end
    %     LogLike = LogLike+log(mean(weights(:,t)))+max(logWeight(:,t));
    %     Loglikelyhood_1(t) = LogLike;

W_T = weights(:,T); W = zeros(N,1); w = zeros(N,1);
[x_T,I] = sort(x_n(:,T));
for i = 1:N
    W(i) = W_T((I(i)));
    w(i) = sum(W(1:i-1));
end
w = w.*(1/w(N));
end
