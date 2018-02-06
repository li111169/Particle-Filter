clc;clear;
rand('seed',1);
N=10;
A=[2,8,2,7,3,5,5,1,4,6];
IndexA=1:N;
W=A./sum(A);
DiedaiNumber=6;
V=[];
AA=A;
WW=W;
for k=1:DiedaiNumber
    outIndex=multinomialR(WW);
    AA=AA(outIndex);
    WW=AA./sum(AA);
    V=[V;AA];
end
figure
subplot(2,1,1);
plot(A','--ro','MarkerFace','g');
subplot(2,1,2);
plot(V(1,:)','--ro','MarkerFace','g');

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