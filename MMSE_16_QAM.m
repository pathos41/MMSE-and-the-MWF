function mmse = MMSE_16_QAM( rou )
%MMSE of BPSK
fun=@(z) (3*exp(-4*rou/5).*sinh(6*sqrt(rou/10).*z)+sinh(2*sqrt(rou/10).*z)).^2.*exp(-z.^2-rou/10)./(exp(-8*rou/10).*cosh(6*sqrt(rou/10).*z)+cosh(2*sqrt(rou/10).*z));
mmse=1-1/(10*sqrt(pi))*integral(fun,-18,18);
end
