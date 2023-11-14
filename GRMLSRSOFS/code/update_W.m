function [W] = update_W(J,M,W,a)
Q = 100;
d = size(W,1);
J1 = a.*eye(d)-J;
for q=1:Q
    Z = 2.*J1*W+2.*M;
    [U,~,V] = svd(Z,"econ");
    W = U*V';
end
end