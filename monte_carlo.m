
clc;clear;
%--------------------------------------------------------------------------
% p=0.5;
% N=1000;
% sum=0;
% for k=1:N
% sum=sum+binornd(1,p);
% P(k)=sum/k;
% end
% figure
% hold on;
% box on;
% plot(1:N,P)
%--------------------------------------------------------------------------

N=5000;m=10000;P=zeros(1,N);
for k=1:N
    P(k)=fun(m);
end
Pave=mean(P);
figure
hold on; box on;
plot(P)
line([0,N],[Pave,Pave],'LineWidth',5,'Color','r');
xlabel('k');
ylabel('estimated probability');

    function p=fun(M)
    frq=0;
    MAX=10;
    for i=1:M
        ball1=unidrnd(MAX);
        ball2=unidrnd(MAX);
        ball3=unidrnd(MAX);
    end
    if((mod(ball1,2)==0)&&(mod(ball2,2)==0)&&(mod(ball3,2)==0))
    frq=frq+1;
    end
    p=frq/M;
    end


