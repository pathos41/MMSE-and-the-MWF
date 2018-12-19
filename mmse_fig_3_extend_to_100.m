snr=-20:20;  %dB scale
rou=10.^(snr/10);
%rou=0:0.1:10;
mmse1=zeros(1,length(rou));
mmse2=zeros(1,length(rou));
mmse3=zeros(1,length(rou));
mmse4=zeros(1,length(rou));
mmse5=zeros(1,length(rou));
for n=1:length(rou)
    fun1=@(x) exp(-x.^2-rou(n))./cosh(2*sqrt(rou(n)).*x);    
    fun2=@(y) exp(-y.^2-rou(n)/2)./cosh(2*sqrt(rou(n)/2).*y);
    
    %a=3*exp(-8*rou(n)/5)*sinh(6*sqrt(rou(n)/5).*z);
    %b=sinh(2*sqrt(rou(n)/5).*z);
    %c=exp(-8*rou(n)/5)*cosh(6*sqrt(rou(n)/5).*z);
    %d=cosh(2*sqrt(rou(n)/5).*z);
    %e=exp(-z.^2-rou(n)/5);
    fun3=@(z) (3*exp(-8*rou(n)/5)*sinh(6*sqrt(rou(n)/5).*z)+sinh(2*sqrt(rou(n)/5).*z)).^2.*exp(-z.^2-rou(n)/5)./(exp(-8*rou(n)/5)*cosh(6*sqrt(rou(n)/5).*z)+cosh(2*sqrt(rou(n)/5).*z));
    
    %f=3*exp(-4*rou(n)/5)*sinh(6*sqrt(rou(n)/10).*u);
    %g=sinh(2*sqrt(rou(n)/10).*u);
    %h=exp(-4*rou(n)/5)*cosh(6*sqrt(rou(n)/10).*u);
    %i=cosh(2*sqrt(rou(n)/10).*u);
    %j=exp(-u.^2-rou(n)/10);
    fun4=@(u) (3*exp(-4*rou(n)/5)*sinh(6*sqrt(rou(n)/10).*u)+sinh(2*sqrt(rou(n)/10).*u)).^2.*exp(-u.^2-rou(n)/10)./(exp(-4*rou(n)/5)*cosh(6*sqrt(rou(n)/10).*u)+cosh(2*sqrt(rou(n)/10).*u));
    
    mmse1(n)=(1/sqrt(pi))*integral(fun1,-inf,inf);   %BPSK
    mmse2(n)=(1/sqrt(pi))*integral(fun2,-inf,inf);   %QPSK
    mmse3(n)=MMSE_4_PAM(rou(n));
    mmse4(n)=MMSE_16_QAM(rou(n));
    %mmse3(n)=1-1/(10*sqrt(pi))*integral(fun3,-inf,inf);  %4-PAM
    %mmse4(n)=1-1/(10*sqrt(pi))*integral(fun4,-inf,inf);  %16-QAM
    mmse5(n)=1/(1+rou(n));             %Gaussian
end

semilogy(snr,mmse1,'b','linewidth',1.25)
hold on
grid on
semilogy(snr,mmse2,'r','linewidth',1.25)
semilogy(snr,mmse3,'k','linewidth',1.25)
semilogy(snr,mmse4,'--','linewidth',1.25)
semilogy(snr,mmse5,'-*')

%ylim([0 1])
xlabel('\rho/dB')
ylabel('MMSE(\rho)')
legend('BPSK','QPSK','4-PAM','16-QAM','Gaussian')