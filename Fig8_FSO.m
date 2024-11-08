clc
clear
syms x y theta
Pt_dBm=25:0.005:35;
Pt=10.^(Pt_dBm./10);
step=0.001;
sigma_X=0.2;
sigma_Y=0.3;
mu_X=0;
mu_Y=-0.1;
A=sqrt(mu_Y.^2+mu_X.^2);
theta_0=pi./2;
gamma_th=1;
A_0=0.8;
h_l=0.01;
w_zeq=1.37;
R_f_low=sqrt(-w_zeq.^2./2.*...
    log(gamma_th./A_0./h_l./Pt));
%% Matlab
CDF_1=[];
theta_0=pi/2;
rho_sym= @(theta)(A.*(cos(theta_0).*cos(theta)./(sigma_X.^2)+...
    sin(theta_0).*sin(theta)./(sigma_Y.^2)));
gamma_0=(cos(theta_0)).^2./(2.*sigma_X.^2)+...
    (sin(theta_0)).^2./(2.*sigma_Y.^2);
gamma_sym= @(theta)((cos(theta)).^2./(2.*sigma_X.^2)+...
    (sin(theta)).^2./(2.*sigma_Y.^2));
for i=1:length(Pt)
    r_a=R_f_low(i);
    part1=@(theta) 1./(2.*gamma_sym(theta)).*...
        (1-exp((-gamma_sym(theta)).*r_a.^2+...
        rho_sym(theta).*r_a));
    part2=@(theta) (rho_sym(theta)).*...
        sqrt(pi)./(4.*(gamma_sym(theta)).^(3/2)).*...
        exp((rho_sym(theta)).^2./...
        (4.*(gamma_sym(theta)))).*...
        (erf((rho_sym(theta))./...
        (2.*sqrt((gamma_sym(theta)))))+....
        erf((2.*(gamma_sym(theta)).*r_a-...
        (rho_sym(theta)))./...
        (2.*sqrt((gamma_sym(theta))))));
    y1=integral(part1,0,2*pi);
    y2=integral(part2,0,2*pi);
    CDF_1(end+1)=1./(2.*pi.*sigma_Y.*sigma_X).*exp...
        (-A.^2.*gamma_0).*(y1+y2);
end
p8=semilogy(Pt_dBm,1-CDF_1,'-*','MarkerIndices',...
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
N=30;
P_LB_N30=[;];
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
        
        P_LB_N30(n,i)=P_LB2T+P_LB1T;
    end
end

P_LB_30=[];
for i=1:length(Pt)
    P_LB_30(i)=0;
    for j=1:1:N-1
        P_LB_30(i)=P_LB_N30(j,i)+P_LB_30(i);
    end
end
p2=semilogy(Pt_dBm,1-P_LB_30,'-s','MarkerIndices',...
    1:600:length(P_LB_30),'color','r');
%%
N=120;
P_LB_N120=[;];
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
        
        P_LB_N120(n,i)=P_LB2T+P_LB1T;
    end
end

P_LB_120=[];
for i=1:length(Pt)
    P_LB_120(i)=0;
    for j=1:1:N-1
        P_LB_120(i)=P_LB_N120(j,i)+P_LB_120(i);
    end
end
p3=semilogy(Pt_dBm,1-P_LB_120,'-','MarkerIndices',...
    1:600:length(P_LB_120),'color','r');




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

N=30;
P_UB_N30=[;];
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
        
        P_UB_N30(n,i)=P_UB2T+P_UB1T;
    end
end
P_UB_30=[];
for i=1:length(Pt)
    P_UB_30(i)=0;
    for j=1:1:N
        P_UB_30(i)=P_UB_N30(j,i)+P_UB_30(i);
    end
end
p5=semilogy(Pt_dBm,1-P_UB_30,'--d','MarkerIndices',...
    1:600:length(P_UB_30),'color','b');

N=120;
P_UB_N3120=[;];
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
        
        P_UB_N3120(n,i)=P_UB2T+P_UB1T;
    end
end
P_UB_120=[];
for i=1:length(Pt)
    P_UB_120(i)=0;
    for j=1:1:N
        P_UB_120(i)=P_UB_N3120(j,i)+P_UB_120(i);
    end
end
p6=semilogy(Pt_dBm,1-P_UB_120,'-.','MarkerIndices',...
    1:600:length(P_UB_120),'color','b');

%%

legend([p8 p1 p2 p3 p4 p5 p6],{'Numerical integration','Upper bound, N=5','Lower bound, N=30',...
    'Upper bound, N=120','Lower bound, N=5',...
    'Upper bound, N=30','Lower bound, N=120',...
    },'Location','southwest','FontSize',13,"FontName","Times New Roman")

axis([25 35 10^(-8) 0.1])



