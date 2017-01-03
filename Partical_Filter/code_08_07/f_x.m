function [x_k_new] = f_x(x_k_old,k)
x_k_new = 0.5*x_k_old+(25*x_k_old)./(1+x_k_old.^2)+8*cos(1.2*k);
end