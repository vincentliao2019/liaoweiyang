%%% 大地电磁二维有限元层状介质matlab程序（TM极化）
%%% 本程序用于矩形单元、双二次插值有限单元法求解二维大地电磁问题
clear
x=[-45 -40 -35 -30 -25 -22.5 -20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10 12.5 15 17.5 20 22.5 25 30 35 40 45];
y=[0 1 2 3 4 5 8 11 14 17 20 23 26 29 32 35 41 50];
A1=zeros(size(x,2)-1,1);
B1=zeros(size(y,2)-1,1);
for i=1:size(x,2)-1
    A1(i)=x(i+1)-x(i);
end
for i=1:size(y,2)-1
    B1(i)=y(i+1)-y(i);
end

Nx=size(A1,1);         %模型的x方向网格数
Ny=size(B1,1);         %模型的y方向网格数
NP=3*Nx*Ny+2*Nx+2*Ny+1;%节点总数
NE=Nx*Ny;              %单元总数
rho=zeros(NE,1);       %各个单元的电阻率
A=zeros(NE,1);
B=zeros(NE,1);

%%%给每个单元的8个节点编号

for IX=1:Nx
    for IY=1:Ny
        N=(IX-1)*Ny+IY;%N为单元编号，(1,N)代表第N个单元的第一个节点
        N1=(IX-1)*(3*Ny+2)+2*IY-1;
        ME(1,N)=N1;
        ME(2,N)=N1+2;
        ME(3,N)=ME(2,N)+3*Ny+2;
        ME(4,N)=N1+3*Ny+2;
        ME(5,N)=N1+1;
        ME(6,N)=ME(2,N)+2*Ny-IY+1;
        ME(7,N)=ME(5,N)+3*Ny+2;
        ME(8,N)=ME(6,N)-1;
        A(N)=A1(IX);
        B(N)=B1(IY);
    end
end

%给模型各个单元赋电阻率的值
for IX=1:10
    for IY=1:Ny
        N=(IX-1)*Ny+IY;            
        rho(N)=10;
    end
end
for IX=11:18
    for IY=1:Ny
        N=(IX-1)*Ny+IY;            
        rho(N)=2;
    end
end
for IX=19:Nx
    for IY=1:Ny
        N=(IX-1)*Ny+IY;            
        rho(N)=1;
    end
end
%%%A和B中存储的值为每个单元的长宽，根据实际区域网格剖分的坐标来确定

f=[0.001 0.005 0.01 0.02 0.025 0.05 0.1 0.125 0.2 0.25 0.5 1 1.25 2 2.5 5 10 12.5 20 25 50 100 500 1000];  %计算的频点
for ff=1:size(f,2)
    K1=zeros(NP);
    K2=zeros(NP);
    K3=zeros(NP);
    P=zeros(NP,1);
    for h=1:NE
        BA=B(h)/A(h)/90;
        AB=A(h)/B(h)/90;
        K=BA*[52 17 23 28  6 -40 -6 -80;17 52 28 23 6 -80 -6 -40;23 28 52 17 -6 -80 6 -40;
            28 23 17 52 -6 -40 6 -80;6 6 -6 -6 48 0 -48 0;-40 -80 -80 -40 0 160 0 80;
            -6 -6 6 6 -48 0 48 0;-80 -40 -40 -80 0 80 0 160]...
            +AB*[52 28 23 17 -80 -6 -40 6;28 52 17 23 -80 6 -40 -6;23 17 52 28 -40 6 -80 -6;
            17 23 28 52 -40 -6 -80 6;-80 -80 -40 -40 160 0 80 0; -6 6 6 -6 0 48 0 -48;
            -40 -40 -80 -80 80 0 160 0;6 -6 -6 6 0 -48 0 48];
        for j=1:8
            NJ=ME(j,h);
            for k=1:8
                NK=ME(k,h);
                K1(NJ,NK)=K1(NJ,NK)+K(j,k)*rho(h);   %rho(h)为单元h的电阻率值
            end
        end
    end
    
   %%%%行程单元K2e %%%%
   mu=(4e-7)*pi;
   w=2*pi*f(ff);
   m=sqrt(-1)*w*mu;
   for h=1:NE
       K=[6 2 3 2 -6 -8 -8 -6;2 6 2 3 -6 -6 -8 -8;3 2 6 2 -8 -6 -6 -8;
           2 3 2 6 -8 -8 -6 -6 ;-6 -6 -8 -8 32 20 16 20;-8 -6 -6 -8 20 32 20 16;
           -8 -8 -6 -6 16 20 32 20;-6 -8 -8 -6 20 16 20 32];
       for j=1:8
           NJ=ME(j,h);
           for k=1:8
               NK=ME(k,h);
               K2(NJ,NK)=K2(NJ,NK)+K(j,k)*m*A(h)*B(h)/180;
           end
       end
   end
   for h1=Ny:Ny:NE
       i=ME(2,h1); j=ME(3,h1);k=ME(4,h1);s=ME(1,h1);
       a=ME(6,h1); b=ME(7,h1);c=ME(8,h1);d=ME(5,h1);
       kn=sqrt(-m*rho(h1));
       mk=kn*A(h1)/30;
       Kjj=4*mk;
       K3(j,j)=K3(j,j)+Kjj;
       K3(i,i)=K3(i,i)+Kjj;
       
       Kji=-1*mk;
        K3(j,i)=K3(j,i)+Kjj;
        K3(i,j)=K3(i,j)+Kjj;
       
       Kai=2*mk;
         K3(a,i)=K3(a,i)+Kai;
         K3(i,a)=K3(i,a)+Kai;
       
         K3(a,j)=K3(a,j)+Kai;
         K3(j,a)=K3(j,a)+Kai;
       
        
       Kaa=16*mk;
         K3(a,a)=K3(a,a)+Kaa;
   end
   
   v=sparse(K1-K2+K3);            %组装总体刚度矩阵
   
   
   
   clear i j k s a b c d mj mk Kii Kjj Kkk Kij Kki Ksi Kaa Kab
   
   %%%%带入上边界值u|AB=1，按照徐世浙院士的《地球物理中的有限元》一书，为了求解方便进行处理
   
   for i=1:Nx+1
       j=1+(i-1)*(3*Ny+2);     %%所有在AB边界上的单元节点
       v(j,j)=v(j,j)*10^10;
       P(j)=v(j,j)*1;
   end
   
   for i=1:Nx
       j=2*Ny+2+(i-1)*(3*Ny+2);   %%所有AB边界上的中点
       v(j,j)=v(j,j)*10^10;
       P(j)=v(j,j);
   end
   
   
   
   %%%求解线性方程组%%%%%
   tol=1e-20;
   [L,U]=lu(v);
   maxit=100;
   an=bicgstab(v,P,tol,maxit,L,U);
   ant(:,ff)=an;
end


%%%%%%%计算视电阻率和阻抗相位%%%%%
vect=1:(3*Ny+2):(NP-2*Ny);
u1=zeros(size(vect,2),size(f,2));
u2=zeros(size(vect,2),size(f,2));
u3=zeros(size(vect,2),size(f,2));
u4=zeros(size(vect,2),size(f,2));
z=zeros(size(vect,2),size(f,2));
pc=zeros(size(vect,2),size(f,2));
phase=zeros(size(vect,2),size(f,2));
l=3*(B(1)/2);
for j=1:size(f,2)
    for i=1:size(vect,2)
        mu=(4e-7)*pi;
        w=2*pi*f(j);
        u1(i,j)=ant(vect(i),j);
        u2(i,j)=ant(vect(i)+1,j);
        u3(i,j)=ant(vect(i)+2,j);
        u4(i,j)=ant(vect(i)+3,j);
        ux(i,j)=(-11*u1(i,j)+18*u2(i,j)-9*u3(i,j)+2*u4(i,j))/(2*l);
        z(i,j)=(rho(1)*ux(i,j))/(ant(vect(i),j));                  %阻抗
        pc(i,j)=(z(i,j)^2)/(w*mu);    %视电阻率
        phase(i,j)=-atan(imag(z(i,j))/real(z(i,j)))*180/pi;      %阻抗相位
    end
end