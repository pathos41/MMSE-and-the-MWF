rou=0:0.1:10;
mmse=zeros(1,length(rou));
mmse1=zeros(1,length(rou));
mmse2=zeros(1,length(rou));
for n=1:length(rou)
    fun=@(y) tanh(2*sqrt(rou(n)/2).*y).*exp(-(y-sqrt(rou(n)/2)).^2)/sqrt(pi);  
    mmse(n)=1-integral(fun,-inf,inf);   %QPSK
    mmse1(n)=1-rou(n);               %Low-power expansion
    mmse2(n)=exp(-rou(n)/2)*(sqrt(pi)-2.1/rou(n))/(sqrt(2)*rou(n));%High-power expansion
end

plot(rou,mmse,'b','linewidth',1.25)
hold on
grid on
plot(rou,mmse1,'-.','linewidth',1.25)
plot(rou,mmse2,'--','linewidth',1.25) 

ylim([0 1])
xlabel('\rho')
ylabel('MMSE(\rho)')
legend('QPSK','low-power expansion','high-power expansion')