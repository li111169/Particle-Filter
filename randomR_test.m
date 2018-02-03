clc;clear;
N=10;
A=[2,8,2,7,3,5,5,1,4,6];
IndexA=1:N;
W=A./sum(A);

OutIndex=randomR(W);

%iteration one
NewA=A(OutIndex);
%iteration two
W=NewA./sum(NewA);
OutIndex=randomR(W);
NewA2=NewA(OutIndex);
%iteration three
W=NewA2./sum(NewA2);
OutIndex=randomR(W);
NewA3=NewA2(OutIndex);

W=NewA3./sum(NewA3);
OutIndex=randomR(W);
NewA4=NewA3(OutIndex);

W=NewA4./sum(NewA4);
OutIndex=randomR(W);
NewA5=NewA4(OutIndex);

W=NewA5./sum(NewA5);
OutIndex=randomR(W);
NewA6=NewA5(OutIndex);

figure
subplot(2,1,1);
plot(A,'--ro','MarkerFace','g');
axis([1,N,1,N])
subplot(2,1,2);
plot(NewA,'--ro','MarkerFace','g');
axis([1,N,1,N]);

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

