
function [x_n,w,x_T] = APF(Y_n,x_initial,N,T,mu,phi,sigma2)
x_n = zeros(N,T); logWeight = zeros(N,T); weights = zeros(N,T);
N_eff = zeros(1,T); x_n_0 = zeros(N,T); Tem = zeros(N,T);
u_k = zeros(N,T);var = zeros(N,T);
Lambda = zeros(N,T); Lambda_new = zeros(N,T);
% LogLike = 0;  Loglikelyhood_1 = zeros(1,T_n

%% set initial value
x_n(:,1) = x_initial;
Tem(:,1) = mu + phi.*(x_n(:,1) - mu);
u_k(:,1) = 0;
var(:,1) = exp(Tem(:,1)+0.5*sigma2);
Lambda(:,1) = -0.5*log(2*pi)-0.5*log(var(:,1))-0.5*((Y_n(1)-u_k(:,1)).^2./var(:,1));
weights(:,1) = (1/N)*ones(N,1);

for t = 2:T
%% assume the weights first

    x_n(:,t) = mu+phi*(x_n(:,t-1)-mu)+sqrt(sigma2)*randn(N,1);
    Tem(:,t) = mu + phi.*(x_n(:,t) - mu);
    u_k(:,t) = 0;
    var(:,t) = exp(u_k(:,t)+0.5*sigma2);
    Lambda(:,t) =-0.5*log(2*pi)-0.5*log(var(:,t))-0.5*((Y_n(t)-u_k(:,t)).^2./var(:,t));
    logWeight(:,t) = Lambda(:,t)+ log(weights(:,t-1));
    weights(:,t) = exp(logWeight(:,t)-max(logWeight(:,t)));
    weights(:,t) = weights(:,t)./sum(weights(:,t));
    
%% Resample
    N_eff(t) = 1/(sum(weights(:,t).^2));
    if N_eff < N
        [~,~,i_j_1] = resample(x_n(:,t),weights(:,t),N);
        for j = 1:N
            x_n_0(j,t-1) = x_n(i_j_1(j),t-1);
            Lambda_new(j,t) = Lambda(i_j_1(j),t);
        end
    else
        x_n_0(:,t-1) = x_n(:,t-1);
        Lambda_new(:,t) = Lambda(:,t);
    end
    
%% Calculate weights again
    x_n(:,t) = mu+phi*(x_n_0(:,t-1)-mu)+sqrt(sigma2)*randn(N,1);
    logWeight(:,t) = -0.5*log(2*pi)-0.5*x_n(:,t)-0.5*(1./exp(x_n(:,t)))*((Y_n(t))^2)+Lambda_new(:,t);
    weights(:,t) = exp(logWeight(:,t)-max(logWeight(:,t)));
    weights(:,t) = weights(:,t)./sum(weights(:,t));
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
