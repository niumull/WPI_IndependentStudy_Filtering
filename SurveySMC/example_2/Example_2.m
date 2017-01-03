clear all
clc
rng(100,'twister')

N = 15; T = 5000; mu = 0.5; phi = 0.985; Y_n = 1;
epsilon_1 = 1; epsilon_2 = 0.5; epsilon_3 = 0.05; 

x_n = zeros(N,T+1); logWeight = zeros(N,T); weights = zeros(N,T);
N_eff = zeros(1,T); x_n_new = zeros(N,T);
% LogLike = 0;  Loglikelyhood_1 = zeros(1,T_n);
for t = 2:T+1
    x_n(:,t) = mu+phi*(x_n(:,t-1)-mu)+randn(N,1);
end
x = x_n(:,2:end);
x_1 = x(1,:);
h_1 = x_1; h_2 = x_1.^2; h_3 = sin(x_1); h_4 = exp(x_1);
w_1_1 = (1/(sqrt(2*pi)*epsilon_1))*exp((-(Y_n-h_1).^2)/(2*epsilon_1^2));
w_1_2 = (1/(sqrt(2*pi)*epsilon_2))*exp((-(Y_n-h_1).^2)/(2*epsilon_2^2));
w_1_3 = (1/(sqrt(2*pi)*epsilon_3))*exp((-(Y_n-h_1).^2)/(2*epsilon_3^2));
w_2_1 = (1/(sqrt(2*pi)*epsilon_1))*exp((-(Y_n-h_2).^2)/(2*epsilon_1^2));
w_2_2 = (1/(sqrt(2*pi)*epsilon_2))*exp((-(Y_n-h_2).^2)/(2*epsilon_2^2));
w_2_3 = (1/(sqrt(2*pi)*epsilon_3))*exp((-(Y_n-h_2).^2)/(2*epsilon_3^2));
w_3_1 = (1/(sqrt(2*pi)*epsilon_1))*exp((-(Y_n-h_3).^2)/(2*epsilon_1^2));
w_3_2 = (1/(sqrt(2*pi)*epsilon_2))*exp((-(Y_n-h_3).^2)/(2*epsilon_2^2));
w_3_3 = (1/(sqrt(2*pi)*epsilon_3))*exp((-(Y_n-h_3).^2)/(2*epsilon_3^2));
w_4_1 = (1/(sqrt(2*pi)*epsilon_1))*exp((-(Y_n-h_4).^2)/(2*epsilon_1^2));
w_4_2 = (1/(sqrt(2*pi)*epsilon_2))*exp((-(Y_n-h_4).^2)/(2*epsilon_2^2));
w_4_3 = (1/(sqrt(2*pi)*epsilon_3))*exp((-(Y_n-h_4).^2)/(2*epsilon_3^2));

[X_1,I] = sort(x_1);
for i = 1:T
    W_1_1(i) = w_1_1(I(i)); W_1_2(i) = w_1_2(I(i)); W_1_3(i) = w_1_3(I(i));
    W_2_1(i) = w_2_1(I(i)); W_2_2(i) = w_2_2(I(i)); W_2_3(i) = w_2_3(I(i));
    W_3_1(i) = w_3_1(I(i)); W_3_2(i) = w_3_2(I(i)); W_3_3(i) = w_3_3(I(i));
    W_4_1(i) = w_4_1(I(i)); W_4_2(i) = w_4_2(I(i)); W_4_3(i) = w_4_3(I(i));
end

figure(1)
subplot(2,2,1);
hold on
plot(X_1,W_1_1,'b-')
plot(X_1,W_1_2,'r--')
plot(X_1,W_1_3,'-.')
legend('sigma_eps = 1','sigma_eps = 0.5','sigma_eps = 0.05','Location','northwest')
xlim([-2.5 2.5])
hold off

subplot(2,2,2);
hold on
plot(X_1,W_2_1,'b-')
plot(X_1,W_2_2,'r--')
plot(X_1,W_2_3,'-.')
legend('sigma_eps = 1','sigma_eps = 0.5','sigma_eps = 0.05','Location','northwest')
xlim([-2.5 2.5])
hold off

subplot(2,2,3);
hold on
plot(X_1,W_3_1,'b-')
plot(X_1,W_3_2,'r--')
plot(X_1,W_3_3,'-.')
legend('sigma_eps = 1','sigma_eps = 0.5','sigma_eps = 0.05','Location','northwest')
xlim([-2.5 2.5])
hold off

subplot(2,2,4);
hold on
plot(X_1,W_4_1,'b-')
plot(X_1,W_4_2,'r--')
plot(X_1,W_4_3,'-.')
legend('sigma_eps = 1','sigma_eps = 0.5','sigma_eps = 0.05','Location','northwest')
xlim([-2.5 2.5])
hold off
