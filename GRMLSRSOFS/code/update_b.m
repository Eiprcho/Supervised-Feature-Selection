function [b] = update_b(W,X,Y,Dv,I1,alpha)
b1 = alpha.*Y'*Dv*I1-W'*X'*Dv*I1;
b = b1./(I1'*Dv*I1);
end

