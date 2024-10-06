function plotpp2(beobachtung,verteilung,parameter1,parameter2)
% PLOTPP Probabilistik-Plot 
%   PLOTPP(BEOBACHTUNG,'VERTEILUNG',PARAMETER1,PARAMETER2)
n=length(beobachtung);
x=sort(beobachtung);

v=(1:n)/n; % Beobachtete Verteilung
if (nargin==2)
    w=cdf(verteilung,x)
elseif (nargin==3)
    w=cdf(verteilung,x,parameter1);
elseif (nargin==4)
    w=cdf(verteilung,x,parameter1, parameter2);
else
    error('PlotPP: Nur bis zu 2 Parameter unterstützt');
end
%mn=min(min(x),min(y));mx=max(max(x),max(y));

% compute critical value for KS statistics
KSLevel=.05;
kolmog=@(x,y)(x<4)*sqrt(2*pi)/x*sum(exp(-(1:2:35).^2*pi^2./(8*x^2)))...
            +(x>=4)*1.0-y;
ks=fzero(@(x)kolmog(x,KSLevel),[0.001,2])/sqrt(n);
%ks=1.2/sqrt(n); %  10%: 1.224, 5%  1.358, 1% 1.628 
plot1=plot(v,w,'.',[0,1],[0,1],'k-',[-ks,1-ks],[0,1],'k:',[ks,1+ks],[0,1],'k:')
set(plot1,'MarkerSize',20)
set(plot1,'LineWidth',2)
set(gca,'FontSize',24)
%plot(v,w,'*','MarkerSize',12,[0,1],[0,1],'k-',[-ks,1-ks],[0,1],'k:','LineWidth',2,[ks,1+ks],[0,1],'k:','LineWidth',2);
axis([0,1,0,1]);
title(['PP Plot']);
xlabel('Observed Probabilities');ylabel('Expected Probabilities');
set(gca,'FontSize',15)