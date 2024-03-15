function [W,obj] = GRMLSRSOFS(X,Y,param)
% Ref: "Sparse orthogonal supervised feature selection with global 
% redundancy minimization, label scaling, and robustness", ins, 2024. 
% Authors: Huming Liao, Hongmei Chen, Yong Mi,Chuan Luo, Shi-Jinn Horng,Tianrui Li.

% Input
% X: d*n data matrix
% Y: n*c label matrix
% param: parameters

%Output
% W: d*c projection matrix
% obj: loss value.

max_ite = 200;
[n,d] = size(X);
c = size(Y,2);
H = eye(n)-(1/n)*ones(n);

alpha = param.alpha;
gama = param.gama;
beta = param.beta;
p = param.p;

I1= ones(n,1);
A = get_A(X,H);
W = zeros(d,c);
b = zeros(c,1);
Dv = eye(n);
D = eye(d);
V = X*W+I1*b'-alpha*Y;
for i=1:max_ite
    %%% calculate the loss
    o1 = 2*trace(V'*Dv*V);
    o2 = 2*gama/p*trace(W'*D*W);
    o3 = beta*trace(W'*A*W);
    obj(i) = o1+o2+o3;
    if i > 1 && abs(obj(i-1)-obj(i)) < 0.000001 %%%reach the convergence condition.
        break;
    end
    
    %%%1.updata b
    b = update_b(W,X,Y,Dv,I1,alpha);

    %%%2.updata W(GPI)
    J = 2.*X'*Dv*X+2*(gama/p).*D+beta.*A;
    M = 2*alpha.*X'*Dv*Y-2.*X'*Dv*I1*b';
    [W] = updata_W(J,M,W,1);

    %%%3.updata D_v,D
    V = X*W+I1*b'-alpha*Y;
    Dv = diag(0.5 * (sqrt(sum(V.^2,2)+eps)).^(-1));
    D = diag(0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2));

end
end