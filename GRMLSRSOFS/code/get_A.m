function [A] = get_A(X,H)
[~,d] = size(X);
for i=1:d
    F(i,:) = H*X(:,i);
end
F = F';
D = diag(1./sqrt(sum(F.^2,1)));
B = (F*D)'*(F*D);
A = B.*B;
end