clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed',1);
T=50;%sampling number
dt=1;%sampling period
Q=1;
R=1;
v=sqrt(R)*randn(T,1);
w=sqrt(Q)*randn(T,1);
numSamples=100;
ResampleStrategy=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=0.1;%initial status

X=zeros(T,1);
Z=zeros(T,1);
X(1,1)=x0;
Z(1,1)=(X(1,1)^2)./20+v(1,1);%Observation intial status

for k=2:T
    %state function
    X(k,1)=0.5*X(k-1,1)+2.5*X(k-1,1)/(1+X(k-1,1)^(2))+8*cos(1.2*k)+w(k-1,1);
    %observation function
    Z(k,1)=(X(k,1).^2)./20+v(k,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle fileter initialize, needs to set up estimation status, particle
%set and weight
Xpf=zeros(numSamples,T);%estimated status
Xparticles=zeros(numSamples,T);%particle set
Zpre_pf=zeros(numSamples,T);%predition for observation
weight=zeros(numSamples,T);%weight initialization

Xpf(:,1)=x0+sqrt(Q)*randn(numSamples,1);
Zpre_pf(:,1)=Xpf(:,1).^2/20;

%update and prediction
%step one
for i=1:numSamples
    QQ=Q;%we do not need keep same variance
    net=sqrt(QQ)*randn;
    Xparticles(i,k)=0.5.*Xpf(i,k-1)+2.5.*Xpf(i,k-1)./(1+Xpf(i,k-1).^2)+8*cos(1.2*k)+net;
end
%step two
for i=1:numSamples
    Zpre_pf(i,k)=Xparticles(i,k)^2/20;
    weight(i,k)=exp(-0.5*R^(-1)*(Z(k,1)-Zpre_pf(i,k))^2);
end
weight(i,k)=weight(:,k)./sum(weight(:,k));
%step three-choose different resampling strategy
if ResampleStrategy==1
    outIndex=randomR(weight(:,k));
elseif ResampleStrategy==2
    outIndex=systematicR(weight(:,k));
elseif ResampleStrategy==3
    outIndex=multinomialR(weight(:,k));
elseif ResampleStrategy==4
    outIndex=residualR(weight(:,k));
end
%step four pick up corresponding particle based on outIndex
Xpf(:,k)=Xparticles(outIndex,k);
end
Xmean_pf=mean(Xpf);
bins=20;
Xmap_pf=zeros(T,1);
for k=1:T
    [p,pos]=hist(Xpf(:,k,1),bins);
    map=find(p==max(p));
    Xmap_pf(k,1)=pos(map(1));
end
for k=1:T
    Xstd_pf(1,k)=std(Xpf(:,k)-X(k,1));
end