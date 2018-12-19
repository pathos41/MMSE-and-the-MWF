%----------------------------Test for fig 9--------------------------------
tol=1e-5;   %Tolerance

c=0;   %Original left side
d=150;  %Original right side

e=0;   %Original left side
f=150;  %Original right side

SNR=-20:20;

Gamma1=10.^(SNR/10);
Gamma2=2.*Gamma1;

p1=zeros(1,length(Gamma1));
p2=zeros(1,length(Gamma1));

for n=1:length(Gamma1)
    gamma1=Gamma1(n);
    gamma2=Gamma2(n);
    a=gamma1*MMSE_QPSK(4*gamma1/3);  %Original left side
    b=gamma2*MMSE_QPSK(4*gamma1/3);  %Original right side
    %max1=-1+ceil((log(b-a)-log(tol))/log(2));   %Number of iterations
    
    %d=Bisection_QPSK(0,100,1e-5,a,gamma2);  %Original right side
    %f=Bisection_QPSK(0,100,1e-5,a,gamma2);  %Original right side
    
    for k=1:20
        eta=(a+b)/2;  %bisection
        
        rou_1a=Bisection_QPSK(c,d,1e-5,a,gamma1);
        rou_2a=Bisection_QPSK(e,f,1e-5,a,gamma2);
        
        rou_1b=Bisection_QPSK(c,d,1e-5,b,gamma1);
        rou_2b=Bisection_QPSK(e,f,1e-5,b,gamma2);
        
        rou_1e=Bisection_QPSK(c,d,1e-5,eta,gamma1);
        rou_2e=Bisection_QPSK(e,f,1e-5,eta,gamma2);
        
        f_a=(1/(2*gamma1))*rou_1a+(1/(2*gamma2))*rou_2a-1;
        f_b=(1/(2*gamma1))*rou_1b+(1/(2*gamma2))*rou_2b-1;
        f_e=(1/(2*gamma1))*rou_1e+(1/(2*gamma2))*rou_2e-1;
        
        if f_e==0   %Find the optimal point
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
            
        else if f_b*f_e>0
                b=eta;
                
            else
                a=eta;
            end
        end
        if b-a<gamma1*MMSE_QPSK(4*gamma1/3)/1e5
            p1(n)=rou_1e/gamma1;
            p2(n)=rou_2e/gamma2;
            break
        end
    end
end



%------------------------High Power Approximation--------------------------

Gamma11=Gamma1;

p11=zeros(1,length(Gamma11));
p21=zeros(1,length(Gamma11));

for n1=1:length(Gamma11)
    gamma11=Gamma11(n1);
    a1=0;  %Original left side
    b1=2;  %Original right side
    max1=-1+ceil((log(b1-a1)-log(tol))/log(2));   %Number of iterations
    
    for k1=1:max1+1
        eta1=(a1+b1)/2;  %bisection
        
        rou_1a1=High_Power(gamma11,a1);
        rou_1b1=High_Power(gamma11,b1);
        rou_1e1=High_Power(gamma11,eta1);
        
        if rou_1e1==0   %Find the optimal point
            p11(n1)=2-eta1;
            p21(n1)=eta1;
            break
            
        else if rou_1b1*rou_1e1>0
                b1=eta1;
                
            else
                a1=eta1;
            end
        end
        if b1-a1<tol
            p11(n1)=2-eta1;
            p21(n1)=eta1;
            break
        end
    end
end

plot(SNR,p1,'linewidth',1.25)
hold on
grid on
plot(SNR,p11,'r','linewidth',1.25)
plot(SNR,p2,'linewidth',1.25)
plot(SNR,p21,'r','linewidth',1.25)

xlabel('h_i^2*P/dB')
ylabel('p_i')
ylim([0 2])  
legend('Exact mercury/waterfilling','High-power Approximation')

text(-5,1.7,'Ch2')
text(-10,1.3,'Ch2')
text(-10,0.7,'Ch1')
text(-5,0.3,'Ch1')