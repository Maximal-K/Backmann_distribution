clc
clear
syms x y theta
Pt_dBm=30:0.005:55;
Pt=10.^(Pt_dBm./10);
step=0.001;

sigma_X=0.1;
sigma_Y=0.15;
mu_X=0.05;
mu_Y=0.1;
A=sqrt(mu_Y^2+mu_X.^2);
theta_0=atan(mu_Y./mu_X);
gamma_th=10;
d_l=5;
D=0.1;
c=0.1;
sigma_o_dl=0.2;
R_f_low=sqrt(-2.*sigma_o_dl.^2.*log(8.*sigma_o_dl.^2./...
   D.^2./Pt.*exp(c.*d_l).*gamma_th));
%%
CDF_1=[];
rho_sym= (A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
    sin(theta_0).*sin(theta)./(sigma_Y.^2)));
gamma_0=(cos(theta_0)).^2./(2.*sigma_X.^2)+...
    (sin(theta_0)).^2./(2.*sigma_Y.^2);
gamma_sym= ((cos(theta)).^2./(2.*sigma_X.^2)+...
    (sin(theta)).^2./(2.*sigma_Y.^2));
for i=1:length(Pt)
    r_a=(R_f_low(i));
    part1=@(theta) 1./(2.*((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2))).*...
        (1-exp(-((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2)).*r_a.^2+...
        (A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
        sin(theta_0).*sin(theta)./(sigma_Y.^2))).*r_a));
    part2=@(theta) (A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
        sin(theta_0).*sin(theta)./(sigma_Y.^2))).*sqrt(pi)./(4.*((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2)).^(3/2)).*...
        exp((A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
        sin(theta_0).*sin(theta)./(sigma_Y.^2))).^2./...
        (4.*((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2)))).*...
        (erf((A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
        sin(theta_0).*sin(theta)./(sigma_Y.^2)))./...
        (2.*sqrt(((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2)))))+....
        erf((2.*((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2)).*r_a-...
        (A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
        sin(theta_0).*sin(theta)./(sigma_Y.^2))))./...
        (2.*sqrt(((cos(theta)).^2./(2.*sigma_X.^2)+...
        (sin(theta)).^2./(2.*sigma_Y.^2))))));
    y1=integral(part1,0,2*pi);
    y2=integral(part2,0,2*pi);
    CDF_1(end+1)=1./(2.*pi.*sigma_Y.*sigma_X).*exp...
        (-A.^2.*gamma_0).*(y1+y2);
end
p8=semilogy(Pt_dBm,1-CDF_1,'-','MarkerIndices',...
    1:600:length(CDF_1),'color','k');
hold on
%%
N=5;
P_LB_N5=[;];
for n=1:1:N-1
    n;
    for i=1:length(Pt)
        
        r_a=R_f_low(i);
        
        P_LB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB_N5(n,i)=P_LB2T+P_LB1T;
    end
end

P_LB_5=[];
for i=1:length(Pt)
    P_LB_5(i)=0;
    for j=1:1:N-1
        P_LB_5(i)=P_LB_N5(j,i)+P_LB_5(i);
    end
end
p1=semilogy(Pt_dBm,1-P_LB_5,'->','MarkerIndices',...
    1:600:length(P_LB_5),'color','r');
hold on
%%
N=20;
P_LB_N20=[;];
for n=1:1:N-1
    n;
    for i=1:length(Pt)
        
        r_a=R_f_low(i);
        
        P_LB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB_N20(n,i)=P_LB2T+P_LB1T;
    end
end

P_LB_20=[];
for i=1:length(Pt)
    P_LB_20(i)=0;
    for j=1:1:N-1
        P_LB_20(i)=P_LB_N20(j,i)+P_LB_20(i);
    end
end
p2=semilogy(Pt_dBm,1-P_LB_20,'-s','MarkerIndices',...
    1:600:length(P_LB_20),'color','r');
%%
N=50;
P_LB_N50=[;];
for n=1:1:N-1
    n;
    for i=1:length(Pt)
        r_a=R_f_low(i);
        P_LB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-(n.*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_LB_N50(n,i)=P_LB2T+P_LB1T;
    end
end

P_LB_50=[];
for i=1:length(Pt)
    P_LB_50(i)=0;
    for j=1:1:N-1
        P_LB_50(i)=P_LB_N50(j,i)+P_LB_50(i);
    end
end
p3=semilogy(Pt_dBm,1-P_LB_50,'-','MarkerIndices',...
    1:600:length(P_LB_50),'color','r');




%%
P_UB_5=[];
N=5;
P_UB_N5=[;];
for n=1:1:N
    n;
    for i=1:length(Pt)
        r_a=R_f_low(i);
        P_UB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB_N5(n,i)=P_UB2T+P_UB1T;
    end
end

for i=1:length(Pt)
    P_UB_5(i)=0;
    for j=1:1:N
        P_UB_5(i)=P_UB_N5(j,i)+P_UB_5(i);
    end
end
p4=semilogy(Pt_dBm,1-P_UB_5,'--<','MarkerIndices',...
    1:600:length(P_UB_5),'color','b');

N=20;
P_UB_N20=[;];
for n=1:1:N
    n;
    for i=1:length(Pt)
        r_a=R_f_low(i);
        P_UB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB_N20(n,i)=P_UB2T+P_UB1T;
    end
end
P_UB_20=[];
for i=1:length(Pt)
    P_UB_20(i)=0;
    for j=1:1:N
        P_UB_20(i)=P_UB_N20(j,i)+P_UB_20(i);
    end
end
p5=semilogy(Pt_dBm,1-P_UB_20,'--d','MarkerIndices',...
    1:600:length(P_UB_20),'color','b');

N=50;
P_UB_N50=[;];
for n=1:1:N
    n;
    for i=1:length(Pt)
        r_a=R_f_low(i);
        P_UB1T=1./4.* ( erf((mu_Y-(n-1).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y-(n).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB2T=1./4.* ( erf((mu_Y+(n).*r_a./N)./sqrt(2)./sigma_Y) -...
            erf((mu_Y+(n-1).*r_a./N)./sqrt(2)./sigma_Y) ) .*...
            ( erf( (mu_X+sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X)-...
            erf( (mu_X-sqrt(r_a.^2-((n-1).*r_a./N).^2))./sqrt(2)./sigma_X) );
        
        P_UB_N50(n,i)=P_UB2T+P_UB1T;
    end
end
P_UB_50=[];
for i=1:length(Pt)
    P_UB_50(i)=0;
    for j=1:1:N
        P_UB_50(i)=P_UB_N50(j,i)+P_UB_50(i);
    end
end
p6=semilogy(Pt_dBm,1-P_UB_50,'-.','MarkerIndices',...
    1:600:length(P_UB_50),'color','b');
xlabel('P_{t}(dBm)',"FontName","Times New Roman")
ylabel('Outage Probability',"FontName","Times New Roman")
legend([p8 p1 p2 p3 p4 p5 p6],{'Numerical integration','Upper bound, N=5','Lower bound, N=20',...
    'Upper bound, N=50','Lower bound, N=5',...
    'Upper bound, N=20','Lower bound, N=50',...
    },'Location','southwest','FontSize',13,"FontName","Times New Roman")
axis([32 55 10^(-5) 0.2])