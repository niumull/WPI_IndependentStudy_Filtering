function [x_k_j,w_k_j,i_j] = resample(x_k,w_k,N)

x_k_j = zeros(1,N);
w_k_j = zeros(1,N);
i_j = 1:N;

c = zeros(1,N);
for i = 2:N
    c(i) = c(i-1)+w_k(i);
end

u = zeros(1,N);
u(1) = rand(1)/N;

n = 1;
for j = 1:N
    u(j) = u(1)+(1/N)*(j-1);
    while u(j)>c(n) && n <= N-1
        n = n+1;
    end
    x_k_j(j) = x_k(n);
    w_k_j(j) = 1/N;
    i_j(j) = n;
end

    
    