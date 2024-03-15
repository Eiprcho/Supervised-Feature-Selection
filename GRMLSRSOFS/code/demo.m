clear,clc,warning('off');
method = "GRMLSRSOFS"; 
load("Yale"); %load data
X = mapminmax(X',0,1)';
Y = full(ind2vec(Y'))';
%set parameters
param.alpha = 1;
param.gama = 1;
param.beta = 1;
param.p = 1;
%train to get W and loss
[W,obj] = GRMLSRSOFS(X,Y,param);
plot(obj(2:end)); %plot the convergence curve