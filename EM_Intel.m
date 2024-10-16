clc
clear all
intel
micro
S=INTEL;
n=length(S);
delta=.02;
par=[.1 .05 .35 .8 .6 .4 .65];

for j=1:20
    fd1=((sqrt(2*pi*par(3)^2*S(1:end-1).^2*delta)).^(-1)).*...
        exp((-((S(2:end)-S(1:end-1)-par(1)*S(1:end-1)).^2))./(2*par(3)^2*S(1:end-1).^2*delta));
    fd2=((sqrt(2*pi*par(4)^2*S(1:end-1).^2*delta)).^(-1)).*...
        exp((-((S(2:end)-S(1:end-1)-par(2)*S(1:end-1)).^2))./(2*par(4)^2*S(1:end-1).^2*delta));
    Pf1(1)=par(7); Pf2(1)=1-par(7);
    for i=2:n-1
        PF1(i)=(Pf1(i-1)*fd1(i))/((Pf1(i-1)*fd1(i))+(Pf2(i-1)*fd2(i)));
        PF2(i)=1-PF1(i);
        Pf1(i)=(par(5)*PF1(i))+((1-par(5))*PF2(i));
        Pf2(i)=1-Pf1(i);
    end
    PS1=zeros(n-1,1)'; PS1(n-1)=PF1(n-1); PS2=zeros(n-1,1)'; PS2(n-1)=PF2(n-1);
    for i=n-2:-1:1
        PS1(i)=((PS1(i+1)*PF1(i)*par(5))/Pf1(i+1))+((PS2(i+1)*PF1(i)*(par(6)))/Pf2(i+1));
        PS2(i)=1-PS1(i);
        if PS1(i)>1
            PS1(i)=1;
        end
        if PS2(i)<0
            PS2(i)=0;
        end
    end
    par1(5)=(sum(PS1(2:end).*(par(5).*PF1(1:end-1)./Pf1(2:end))))/(sum(PS1(1:end-1)));
    par1(6)=(sum(PS2(2:end).*(par(6).*PF1(1:end-1)./Pf2(2:end))))/(sum(PS1(1:end-1)));
    par1(7)=PS1(1);
    
    S1=(S(2:end)-S(1:end-1))./S(1:end-1);
    par1(1)=sum(PS1.*S1')/(sum(PS1));
    par1(2)=sum(PS2.*S1')/(sum(PS2));
    S2=((S(2:end)-S(1:end-1)-par1(1)*S(1:end-1)).^2)./(S(1:end-1).^2*delta);
    S3=((S(2:end)-S(1:end-1)-par1(2)*S(1:end-1)).^2)./(S(1:end-1).^2*delta);
    par1(3)=sqrt((sum(PS1.*S2'))/(sum(PS1)));
    par1(4)=sqrt((sum(PS2.*S2'))/(sum(PS2)));
    par=par1;
    r1(j)=par(1); r2(j)=par(2);
    sigma1(j)=par(3); sigma2(j)=par(4);
    p11(j)=par(5); p12(j)=1-par(5);
    p21(j)=par(6); p22(j)=1-par(6);
end
figure(1)
clf
par(1:6)'
subplot(4,2,1)
plot(r1)
subplot(4,2,2)
plot(r2)
subplot(4,2,3)
plot(sigma1)
subplot(4,2,4)
plot(sigma2)
subplot(4,2,5)
plot(p11)
subplot(4,2,6)
plot(p12)
subplot(4,2,7)
plot(p21)
subplot(4,2,8)
plot(p22)








%%%%%%%%%%%%%%






for i=1:n-1
    if PS1(i)>.5
        r(i)=par(1);
        sigma(i)=par(3);
    else
        r(i)=par(2);
        sigma(i)=par(4);
    end
end
RSGBM(1)=S(1);
RSGBM1(1)=S(1);
for i=2:n
    RSGBM(i)=RSGBM(i-1)+RSGBM(i-1)*r(i-1)*delta+sigma(i-1)*RSGBM(i-1)*sqrt(delta)*randn;
    RSGBM1(i)=S(i-1)+S(i-1)*r(i-1)*delta+sigma(i-1)*S(i-1)*sqrt(delta)*randn;
end

figure(2)
clf
hold on
plot(S(1),'b','linewidth',2)
plot(S(2),'r-.','linewidth',2)
legend('regime 1', 'regime 2')
for i=1:n-1
    hold on
    if sigma(i)==par(3)
        plot([i i+1],[RSGBM1(i) RSGBM1(i+1)],'Color','b','linewidth',1.5)
    else
       plot([i i+1],[RSGBM1(i) RSGBM1(i+1)],'r-.','linewidth',1)
    end
end





