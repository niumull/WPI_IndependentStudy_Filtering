clear all
clc 

rng(100,'twister')

N = 15; N_1 = 15; N_2 = 250; N_3 = 1000; N_4 = 10000; T = 50; t_1 = 5; 
t_2 = 10; t_3 = 30; t_4 = 50;  mu = 0.5; phi = 0.985; sigma2 = 0.04;

% True value
x_n = zeros(1,T+100);
for t=1:T+100
    if(t==1)
        x_n(t) = randn(1,1);
    else
        x_n(t) = mu + phi*(x_n(t-1) - mu) + sqrt(sigma2)*randn(1);
    end
end
x_n_true = x_n(101:end);
Y_n = exp(x_n_true./2).*randn(1,T);
% x_initial = mu+sqrt(sigma2/(1-phi^2))*randn(N,1);
% x_initial = randn(1)*ones(N,1);
x_initial = x_n_true(1);

[x_n_1,~,~] = SIR(Y_n,x_initial,N,t_1,mu,phi,sigma2);
[x_n_2,~,~] = SIR(Y_n,x_initial,N,t_2,mu,phi,sigma2);
[x_n_3,~,~] = SIR(Y_n,x_initial,N,t_3,mu,phi,sigma2);
[x_n_4,~,~] = SIR(Y_n,x_initial,N,t_4,mu,phi,sigma2);

% x_initial_1 = mu+sqrt(sigma2/(1-phi^2))*randn(N_1,1);
% x_initial_2 = mu+sqrt(sigma2/(1-phi^2))*randn(N_2,1);
% x_initial_3 = mu+sqrt(sigma2/(1-phi^2))*randn(N_3,1);
% x_initial_4 = mu+sqrt(sigma2/(1-phi^2))*randn(N_4,1);
x_initial_1 = randn(1)*ones(N_1,1);
x_initial_2 = randn(1,1)*ones(N_2,1);
x_initial_3 = randn(1,1)*ones(N_3,1);
x_initial_4 = randn(1,1)*ones(N_4,1);

[~, w_1, x_T_1] = SIR(Y_n,x_initial_1,N_1,T,mu,phi,sigma2);
[~, w_2, x_T_2] = SIR(Y_n,x_initial_2,N_2,T,mu,phi,sigma2);
[~, w_3, x_T_3] = SIR(Y_n,x_initial_3,N_3,T,mu,phi,sigma2);
[~, w_4, x_T_4] = SIR(Y_n,x_initial_4,N_4,T,mu,phi,sigma2);


figure(1)
subplot(2,2,1);
plot(x_n_true,'b-')
hold on
plot(x_n_1','r-+')
legend('true log-volatility','Location','northwest')
xlim([0 50])
ylim([-3 3])
hold off

subplot(2,2,2);
plot(x_n_true,'b-')
hold on
plot(x_n_2','r-+')
legend('true log-volatility','Location','northwest')
xlim([0 50])
ylim([-3 3])
hold off

subplot(2,2,3);
plot(x_n_true,'b-')
hold on
plot(x_n_3','r-+')
legend('true log-volatility','Location','northwest')
xlim([0 50])
ylim([-3 3])
hold off

subplot(2,2,4);
plot(x_n_true,'b-')
hold on
plot(x_n_4','r-+')
legend('true log-volatility','Location','northwest')
xlim([0 50])
ylim([-3 3])
hold off

figure(2)
subplot(2,2,1);
stairs(x_T_1,w_1,'r-')
legend('(i)N = 15','Location','northwest')

subplot(2,2,2);
stairs(x_T_2,w_2,'r-')
legend('(ii)N = 250','Location','northwest')

subplot(2,2,3);
stairs(x_T_3,w_3,'r-')
legend('(iii)N = 1000','Location','northwest')

subplot(2,2,4);
stairs(x_T_4,w_4,'r-')
legend('(iv)N = 10000','Location','northwest')

