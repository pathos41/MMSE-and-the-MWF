function mmse = MMSE_4_PAM( rou )
%MMSE of BPSK
fun=@(z) (3*exp(-8*rou/5).*sinh(6*sqrt(rou/5).*z)+sinh(2*sqrt(rou/5).*z)).^2.*exp(-z.^2-rou/5)./(exp(-8*rou/5).*cosh(6*sqrt(rou/5).*z)+cosh(2*sqrt(rou/5).*z));
mmse=1-1/(10*sqrt(pi))*integral(fun,-18,18);
end
