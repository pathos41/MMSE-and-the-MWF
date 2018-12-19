gamma1=2;
gamma2=2;
gamma3=2;
p=0:0.01:3;

mmse1=zeros(1,length(p));
mmse2=zeros(1,length(p));
mmse3=zeros(1,length(p));

a=zeros(1,length(p));
b=zeros(1,length(p));
c=zeros(1,length(p));
eta=zeros(1,length(p));

rou1=zeros(1,length(p));
rou2=zeros(1,length(p));

syms x y

for n=1:length(p)
    rou1(n)=p(n)*gamma1;
    rou2(n)=p(n)*gamma2;
    
    fun1=@(x) tanh(2*sqrt(rou1(n)).*x).*exp(-(x-sqrt(rou1(n))).^2)/sqrt(pi);   %BPSK
    fun2=@(y) tanh(2*sqrt(rou2(n)/2).*y).*exp(-(y-sqrt(rou2(n)/2)).^2)/sqrt(pi);   %QPSK
    
    mmse1(n)=1-integral(fun1,-inf,inf);
    mmse2(n)=1-integral(fun2,-inf,inf);
    mmse3(n)=1/(1+gamma3*p(n));  %Gaussian
    
    a(n)=gamma1*mmse1(n);
    b(n)=gamma2*mmse2(n);
    c(n)=gamma3*mmse3(n);
    eta(n)=0.48;
end

plot(p,a,'b','linewidth',1.25)
hold on
grid on
plot(p,b,':','linewidth',1.25)
plot(p,c,'-.','linewidth',1.25)
plot(p,eta,'--','linewidth',1.25)

xlabel('p_i')
ylabel('\gamma_iMMSE(p_i\gamma_i)')
legend('BPSK','QPSK','Gaussian','\eta')