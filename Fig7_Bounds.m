clc;clear;
close
r=0:0.001:1;
mu_X=1/2;
sigma_X=1;
mu_Y=1;
sigma_Y=1/2;
%%
%exact CDF(对（4）式符号积分)
k=@(x,y)1./2./pi./sigma_X./sigma_Y.*exp(-(x-mu_X).^2./2./sigma_X./sigma_X-(y-mu_Y).^2./2./sigma_Y./sigma_Y);
F_R4=[];
for i=1:1001
min=@(x)-sqrt(r(i).^2-x.^2);
max=@(x)sqrt(r(i).^2-x.^2);
F_R4(i)=integral2(k,-r(i),r(i),min,max);
end
%%
%下界
FL1=lb(r,2,mu_X,sigma_X,mu_Y,sigma_Y);
FL2=lb(r,5,mu_X,sigma_X,mu_Y,sigma_Y);
FL3=lb(r,10,mu_X,sigma_X,mu_Y,sigma_Y);
FL4=lb(r,20,mu_X,sigma_X,mu_Y,sigma_Y);

semilogy(r,F_R4,'k',r,FL1,'y',r,FL2,'g',r,FL3,'b',r,FL4,'r');

hold on
grid on

%上界
FU1=ub(r,2,mu_X,sigma_X,mu_Y,sigma_Y);
FU2=ub(r,5,mu_X,sigma_X,mu_Y,sigma_Y);
FU3=ub(r,10,mu_X,sigma_X,mu_Y,sigma_Y);
FU4=ub(r,20,mu_X,sigma_X,mu_Y,sigma_Y);

semilogy(r,FU1,'y',r,FU2,'g',r,FU3,'b',r,FU4,'r');
%%
ylim([4.46e-6 0.5])
xlim([-0.02 0.95])
title('Fig.7. CDF bounds of Beckmann distribution and the exact CDF using numerical integral estimation.The step length for the numerical estimation is 0.001.')
xlabel('r')
ylabel('F_R(r)')
legend('F_R(r)','Bounds on F_R(r),N=2','Bounds on F_R(r),N=5','Bounds on F_R(r),N=10','Bounds on F_R(r),N=20','Location','northwest')

axes('Position',[0.545,0.18,0.3,0.4]); % 生成子图   左右  上下 宽窄 
semilogy(r,F_R4,'k',r,FL1,'y',r,FL2,'g',r,FL3,'b',r,FL4,'r');
hold on
grid on
semilogy(r,FU1,'y',r,FU2,'g',r,FU3,'b',r,FU4,'r');
xlim([0.615 0.6265])
ylim([0.02 0.09])
xlabel('r')
ylabel('F_R(r)')

function F=lb(r,N,u1,s1,u2,s2)
F=0;
for n=1:1:N-1
    F=F+1./4.*(erf((u2-(n-1).*r./N)./sqrt(2)./s2)-erf((u2-n.*r./N)./sqrt(2)./s2)).*(erf((u1+sqrt(r.^2-(n.*r./N).^2))./sqrt(2)./s1)-erf((u1-sqrt(r.^2-(n.*r./N).^2))./sqrt(2)./s1))...
        +1./4.*(erf((u2+n.*r./N)./sqrt(2)./s2)-erf((u2+(n-1).*r./N)./sqrt(2)./s2)).*(erf((u1+sqrt(r.^2-(n.*r./N).^2))./sqrt(2)./s1)-erf((u1-sqrt(r.^2-(n.*r./N).^2))./sqrt(2)./s1));
end
end


function F=ub(r,N,u1,s1,u2,s2)
F=0;
for n=1:1:N
    F=F+1./4.*(erf((u2-(n-1).*r./N)./sqrt(2)./s2)-erf((u2-n.*r./N)./sqrt(2)./s2)).*(erf((u1+sqrt(r.^2-((n-1).*r./N).^2))./sqrt(2)./s1)-erf((u1-sqrt(r.^2-((n-1).*r./N).^2))./sqrt(2)./s1))...
        +1./4.*(erf((u2+n.*r./N)./sqrt(2)./s2)-erf((u2+(n-1).*r./N)./sqrt(2)./s2)).*(erf((u1+sqrt(r.^2-((n-1).*r./N).^2))./sqrt(2)./s1)-erf((u1-sqrt(r.^2-((n-1).*r./N).^2))./sqrt(2)./s1));
end
end