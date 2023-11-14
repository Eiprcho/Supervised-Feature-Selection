function [W] = GRMLSRSOFS(X,Y,param)
max_ite = 30;
[n,d] = size(X);
c = size(Y,2);
H = eye(n)-(1/n)*ones(n);

alpha = param.alpha;
gama = param.gama;
beta = param.beta;
p = param.p;
a = param.a;

I1= ones(n,1);
A = get_A(X,H);
W = zeros(d,c);
b = zeros(c,1);
Dv = eye(n);
D = eye(d);
V = X*W+I1*b'-alpha*Y;
for i=1:max_ite

    %%%1.update b
    b = update_b(W,X,Y,Dv,I1,alpha);

    %%%2.update W(GPI)
    J = 2.*X'*Dv*X+2*(gama/p).*D+beta.*A;
    M = 2*alpha.*X'*Dv*Y-2.*X'*Dv*I1*b';
    [W] = update_W(J,M,W,a);

    %%%3.update D_v,D
    V = X*W+I1*b'-alpha*Y;
    Dv = diag(0.5 * (sqrt(sum(V.^2,2)+eps)).^(-1));
    D = diag(0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2));

end
end