close all
clear all
clc

r12=0.5;
r13=0.5;
r23=0.5;
muX=[1 1 1 1]
SigX=[1 r12 r13 0; r12 1 r23 0; r13 r23 1 0; 0 0 0 1];
b=[0 0]';
A=[4 -2 1 0; 2 5 -1 0];
n=8192;k=4;
%u=net(sobolset(3),n+1); u(1,:)=[];
u=sobolpoints(n,k);
%u=rand(n,k);
X=norminv(u)*sqrtm(SigX)+muX; % chol()
Y=(A*X'+b)';

%
nUsed=n;
%nUsed=n;

Si_ed_kern=mmd2si(X,Y,32,'ed')*2
% out-of-memory for whole sample
Si_ed_direct=edistmim(X(1:(nUsed),:),Y(1:(nUsed),:),32)*2

%%
clear u X Y
n=10000;
R=400;

Sis=zeros(R,k);
SSig=sqrtm(SigX); % chol()
if isempty(gcp('nocreate')), parpool('threads');end
parfor r=1:R
%    if(mod(r,10)==0), disp(r); end
  x=norminv(rand(n,k))*SSig+muX; % chol()
  y=(A*x'+b)';
  res=mmd2si(x,y,22,'ed')*2;
  Sis(r,:)=res;
end
figure
for j=k:k
    subplot(1,2,1)
    histogram(Sis(:,j),50);
    xlabel(['x_{' num2str(j) '}']);
    title('Histogram')
    set(gca,'FontSize',15)
    subplot(1,2,2);
    plotpp2(Sis(:,j),'norm',mean(Sis(:,j)),std(Sis(:,j)));
end
