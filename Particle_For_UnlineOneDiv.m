clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('seed',1);
T=50;%sampling number
dt=1;%sampling period
Q=1;
R=1;
v=sqrt(R)*randn(T,1);
w=sqrt(Q)*randn(T,1);
numSamples=100;
ResampleStrategy=1;
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle fileter initialize, needs to set up estimation status, particle
%set and weight
Xpf=zeros(numSamples,T);%estimated status
Xparticles=zeros(numSamples,T);%particle set
Zpre_pf=zeros(numSamples,T);%predition for observation
weight=zeros(numSamples,T);%weight initialization

Xpf(:,1)=x0+sqrt(Q)*randn(numSamples,1);
Zpre_pf(:,1)=Xpf(:,1).^2/20;
for k=2:T
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
weight(:,k)=weight(:,k)./sum(weight(:,k));

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
end

%step four pick up corresponding particle based on outIndex
Xpf(:,k)=Xparticles(outIndex,k);
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
figure(1);clf;
subplot(221);
plot(v);
xlabel('time');
ylabel('measured noise','fontsize',15);
subplot(222);
plot(w);
xlabel('time');
ylabel('process noise','fontsize',15)
subplot(223);
plot(X);
xlabel('time','fontsize',15);
ylabel('real status x','fontsize',15);
subplot(224);
plot(Z);
xlabel('time','fontsize',15);
ylabel('observation Z','fontsize',15);
figure(2);clf;
k=1:dt:T;
plot(k,X,'b',k,Xmean_pf,'r',k,Xmap_pf,'g');
legend('real system status','posterior estimation','maximum posterior estimation');
xlabel('time','fontsize',15);
ylabel('status estimation','fontsize',15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outIndex=randomR(weight)
L=length(weight);
outIndex=zeros(1,L);
u=unifrnd(0,1,1,L);
u=sort(u);
cdf=cumsum(weight);
i=1;
for j=1:L
    while(i<=L)&&(u(i)<=cdf(j))
        outIndex(i)=j;
        i=i+1;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outIndex=systematicR(weight)
% [Row,N]=size(weight);
N=length(weight);
N_children=zeros(1,N);
label=zeros(1,N);
label=1:1:N;
s=1/N;
auxw=0;
auxl=0;
li=0;
T=s*rand(1);
j=1;
Q=0;
i=0;
u=rand(1,N);
while(T<1)
    if(Q>T)
        T=T+s;
        N_children(1,li)=N_children(1,li)+1;
    else
        i=fix((N-j+1)*u(1,j))+j;
        auxw=weight(1,i);
        li=label(1,i);
        Q=Q+auxw;
        weight(1,i)=weight(1,j);
        label(1,i)=label(1,j);
        j=j+1;
    end
end
index=1;
for i=1:N
    if(N_children(1,i)>0)
        for j=index:index+N_children(1,i)-1
            outIndex(j)=i;
        end
end
index=index+N_children(1,i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outIndex=multinomialR(weight)
Col=length(weight);
N_babies=zeros(1,Col);
cdf=cumsum(weight);
u=rand(1,Col);
uu=u.^(1./(Col:-1:1));
ArrayTemp=cumprod(uu);
u=fliplr(ArrayTemp);
j=1;
for i=1:Col
    while(u(i)>cdf(j))
        j=j+1;%position shift
    end
    N_babies(j)=N_babies(j)+1;
end
index=1;
for i=1:Col
    if(N_babies(i)>0)
        for j=index:index+N_babies(i)-1
            outIndex(j)=i;
        end
    end
    index=index+N_babies(i);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outIndex=residualR(weight)
N=length(weight);
N_babies=zeros(1,N);
q_res=N.*weight;
N_babies=fix(q_res);
N_res=N-sum(N_babies);
if(N_res~=0)
    q_res=(q_res-N_babies)/N_res;
    cumDist=cumsum(q_res);
    u=fliplr(cumprod(rand(1,N_res).^(1./(N_res:-1:1))));
    j=1;
    for i=1:N_res
        while(u(1,i)>cumDist(1,j))
            j=j+1;
        end
        N_babies(1,j)=N_babies(1,j)+1;
    end
end
index=1;
for i=1:N
    if(N_babies(1,i)>0)
        for j=index:index+N_babies(1,i)-1
            outIndex(j)=i;
        end
    
    end
    index=index+N_babies(1,i);
end
end