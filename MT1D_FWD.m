function [rho_a,phase]=MT1D_FWD(rho,h)
%�������rhoΪ���������ֵ��hΪ���ֵ
%�������Ǵ�ص�ų����״���ʽ�����ļ��㣬������Ϊ��ά����Ԫ����״���ʵĽ���Աȳ���
mu=(4e-7)*pi;
T=logspace(-3,4,40);
k=zeros(size(rho,2),size(T,2));
for N=1:size(rho,2)
    k(N,:)=sqrt(-i*2*pi*mu./(T.*rho(N)));
end
m=size(rho,2);
Z=-(i*mu*2*pi)./(T.*k(m,:));
for n= m-1:-1:1
    A=-(i*2*pi*mu)./(T.*k(n,:));
    B=exp(-2*k(n,:)*h(n));
    Z=A.*(A.*(1-B)+Z.*(1+B))./(A.*(1+B)+Z.*(1-B));
end
rho_a=(T./(mu*2*pi)).*(abs(Z).^2);
phase=-atan(imag(Z)./real(Z)).*180/pi;



%������һ��ģ����֤�ľ��������Ը����Լ���Ҫ���и�д
% T=logspace(-3,3,40);
% [rho_a,phase]=MT1D_FWD([100,1000,500,10],[1000,2000,3400,130000]);
% subplot(2,1,1)
% loglog(T,rho_a,'-*')
% xlabel('T(s)')
% ylabel('\rho_a(\Omega\cdotm)')
% subplot(2,1,2)
% semilogx(T,phase,'-*')
% xlabel('T(s)')
% ylabel('phase(\circ)')