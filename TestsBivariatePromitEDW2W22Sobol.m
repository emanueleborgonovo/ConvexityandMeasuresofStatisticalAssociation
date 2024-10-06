close all
clear all
clc
if isempty(gcp('nocreate'))==1
parpool("threads")
end
tic
 n=2^15; % 8192;
 k=4;
 R=1;
 res=[];
 for r=1:R
  fprintf('Run %d of %d\n',r,R);
  %u=rand(n,k);
  u=sobolpoints(n,k);
  x=[1,1,1,0]+norminv(u)*sqrtm([1 .5 .5 0; .5 1 .5 0; .5 .5 1 0; 0 0 0 1]);
  y=x*[4 -2 1 0;  2 5 -1 0]';%[4 -2 1 0]'; % [ 4 -2 1 0]?  % ;
  Ms=[5:5:100];
  Vy=sum(var(y));
  for M=Ms
    M
    [Ei,Si]=mmd2si(x,y,M,'ed');
    [W2,W22sep]=bwsi(x,y,M);
    res(end+1,:)=[r,M,Si',Ei,W2,mean(W22sep)/(2*Vy)];
  end
 end
 %%
 for i=1:R
  subplot(2,2,1);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+(1:k)),'o-')
title('Sobol');
  subplot(2,2,2);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+k+(1:k)),'o-')
title('Energy Distance');
  subplot(2,2,3);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+2*k+(1:k)),'o-')
  title('Wasserstein2');
  subplot(2,2,4);
  set(gca,'ColorOrderIndex',1);
  semilogx(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+3*k+(1:k)),'o-')
  title('Wasserstein22');
 end
 toc
%% please change
edist_ana=[3.47;	3.63;	0.94;0];
w22_ana=[6.47;	6.52;	2.86;0];
w2_ana=[6.47; 6.52; 2.86;0];
 figure
  for i=1:R
  subplot(1,2,1);
  set(gca,'ColorOrderIndex',1);%semilogx
  plot(Ms,2*res((i-1)*length(Ms)+(1:length(Ms)),2+k+(1:k)),'o-','LineWidth',2); % energy dist
  title('Estimates of \xi^{ED}(Y,X_i)','FontSize',24);xlabel('Partition Size','FontSize',24)
  ylabel('\xi^{ED}(Y,X_i)','FontSize',24);
  hold on
  a=axis;
  plot([a(1) a(2)],edist_ana*[1,1],'k:','LineWidth',2);%semilogx
  ax=gca;
  ax.XLabel.FontSize=20; 
  legend('$\widehat{\xi}^{ED}(Y,X_1)$','$\widehat{\xi}^{ED}(Y,X_2)$','$\widehat{\xi}^{ED}(Y,X_3)$','$\widehat{\xi}^{ED}(Y,X_4)$','FontSize',24,'Interpreter','Latex')
  hold off
 % subplot(2,2,3);
 % loglog(Ms,sum(abs( res((i-1)*length(Ms)+(1:length(Ms)),2+k+(1:k)) - edist_ana'),2),'--');
  subplot(1,2,2);
  set(gca,'ColorOrderIndex',1);
  plot(Ms,res((i-1)*length(Ms)+(1:length(Ms)),2+2*k+(1:k)),'o-','LineWidth',2); % wasserstein2 (sqrt); semilogx
  title('Estimates of \xi^{W2}(Y,X_i)','FontSize',24);xlabel('Partition Size','FontSize',24)
  ylabel('\xi^{W2}(Y,X_i)','FontSize',24);
    hold on
  a=axis;
  plot([a(1) a(2)],w22_ana*[1,1],'k:','LineWidth',2);%semilogx
   legend('$\widehat{\xi}^{W2}(Y,X_1)$','$\widehat{\xi}^{W2}(Y,X_2)$','$\widehat{\xi}^{W2}(Y,X_3)$','$\widehat{\xi}^{W2}(Y,X_4)$','FontSize',24,'Interpreter','Latex')

  ax=gca;
  ax.XLabel.FontSize=20; 
  hold off
  end

  %save('PromitResults.mat','res','-mat')
  %save('PromitResults.txt','res','-ascii')
%
 %%
delete(gcp)