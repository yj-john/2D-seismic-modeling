%2D acoustic modeling using staggered grid
%����������������ģ��
%ʮ�׽������񾫶�
clc;
clear;
tic;
X=1000,Y=1000;%ģ�ͼ��������С
DX=10,DY=10;%�������С
Npad=20;
NX=X/DX,NY=Y/DY;%����ʵ����
nx=NX+10+2*Npad,ny=NY+10+2*Npad;%���������Ĵ�С
SX=500,SY=500;%��Դλ��
sx=SX/DX+5+Npad,sy=SY/DY+5+Npad;%%��Դ����λ��
Vp=1500*ones(nx,ny);
e=2*ones(nx,ny);
nt=0.6%����ģ��ʱ��
dt=0.001%����ʱ����
ft=40%��Ƶ��С
%�׿��Ӳ�
ns=100;
it=(1:100)*dt;
wave=(1-2*(pi*ft.*(it-1/ft).^2)).*exp(-(pi*ft.*(it-1/ft)).^2);
figure(1); 
plot(1:ns,wave);xlabel('�������');
%������������
%����������ϵ��
fdcoeff=[1.1211  -0.089722 0.013843 -0.0017657 0.00011868]
isnap=1;
%������ʼ��
vx1=zeros(nx,ny);%ǰһ��ʱ��(k)��Vx����
vx2=zeros(nx,ny);%��һ��ʱ��(k-1)��Vx����
vy1=zeros(nx,ny);%ǰһ��ʱ��(k)��Vy����
vy2=zeros(nx,ny);%��һ��ʱ��(k-1)��Vy����
p=zeros(nx,ny);
p1=zeros(nx,ny);
px=zeros(nx,ny);%ǰһ��ʱ��(k+1/2)�ķ���
py=zeros(nx,ny);
px1=zeros(nx,ny);%��һ��ʱ��(k-1/2)�ķ���
py1=zeros(nx,ny);
dx=zeros(nx,ny);
dy=zeros(nx,ny);
%�������ճ�ʼ��
temp=zeros(NX,NY);
%���ձ߽�ϵ����ʼ��
 %��ȫƥ���  ��˥��ģ�͹�ʽd1(x)=3*Vmax/(2*npad)*log(1/R)(x/npad)^2
 %d2(z)=3*Vmax/(2*npad)*log(1/R)(z/npad)^2
 % left //���MPL�߽�˥��
     for i=6:Npad+5
         for j=Npad+6:(ny-(Npad+5))
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-i)/Npad]^2;
             dy(i,j)=0;
         end
     end
% right
    for i=(nx-Npad-4):(nx-5)
         for j=Npad+6:(ny-(Npad+5))
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(i-(nx-5-Npad))/Npad]^2;
             dy(i,j)=0;
         end
    end
% up ���ϱ߽�˥��
   for i=Npad+6:(nx-(Npad+5))
       for j=6:Npad+5
           dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-j)/Npad]^2;
           dx(i,j)=0;
       end
   end
%down ���±߽�˥��
   for i=Npad+6:(nx-(Npad+5))
      for j=(ny-Npad-4):(ny-5)
            dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(j-(ny-5-Npad))/Npad]^2;
            dx(i,j)=0;
      end
   end
%�ĸ��߽�
     for  i=6:Npad+5
         for j=6:Npad+5
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-i)/Npad]^2;
             dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-j)/Npad]^2;
         end
     end
     for  i=6:Npad+5
         for j=(ny-Npad-4):(ny-5)
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-i)/Npad]^2;
             dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(j-(ny-5-Npad))/Npad]^2;
         end
     end
     for  i=(nx-Npad-4):(nx-5)
         for j=6:Npad+5
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(i-(nx-5-Npad))/Npad]^2;
             dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-j)/Npad]^2;
         end
     end
     for  i=(nx-Npad-4):(nx-5)
         for j=(ny-Npad-4):(ny-5)
             dx(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(i-(nx-5-Npad))/Npad]^2;
             dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(j-(ny-5-Npad))/Npad]^2;
         end
     end
%��������
for n=1:nt/dt
dvx=zeros(nx,ny);%Vx��������
dvy=zeros(nx,ny);%Vy��������
dpx=zeros(nx,ny);%dp��x������
dpy=zeros(nx,ny);%dp��y������

    for i=6:(nx-5)
       for j=1:ny
          for k=1:5
            dpx(i,j)=fdcoeff(k)*(p1(i+k,j)-p1(i+1-k,j))+dpx(i,j);
          end
          vx1(i,j)=(1-dt*dx(i,j))*vx2(i,j)-(dt/(DX*e(i,j)))*dpx(i,j);
       end
    end
   for i=1:nx
        for j=6:(ny-5) 
           for k=1:5
            dpy(i,j)=fdcoeff(k)*(p1(i,j+k)-p1(i,j+1-k))+dpy(i,j);
           end
            vy1(i,j)=(1-dt*dy(i,j))*vy2(i,j)-(dt/(DY*e(i,j)))*dpy(i,j);
        end
   end 
    for i=6:(nx-5)
        for j=1:ny
            for k=1:5
            dvx(i,j)=fdcoeff(k)*(vx1(i-1+k,j)-vx1(i-k,j))+dvx(i,j);
            end
            px(i,j)=(1-dt*dx(i,j))*px1(i,j)-dt*e(i,j)*Vp(i,j)^2*(1.0/DX)*(dvx(i,j));
        end
    end
   for i=1:nx
        for j=6:(ny-5)
            for k=1:5
            dvy(i,j)=fdcoeff(k)*(vy1(i,j-1+k)-vy1(i,j-k))+dvy(i,j);
            end
            py(i,j)=(1-dt*dy(i,j))*py1(i,j)-dt*e(i,j)*Vp(i,j)^2*(1.0/DY)*(dvy(i,j));
        end
   end
    for i=1:nx
        for j=1:ny
            p(i,j)=(px(i,j)+py(i,j));
        end
    end
 %��Դ��������
    if n<ns
        p(sx,sy)=wave(n);
    end
      %��������
    p1=p;
    px1=px;
    py1=py;
   vx2=vx1;
   vy2=vy1;
    %�����������
    for i=1:NX
        for j=1:NY
            temp(i,j)=p(i+Npad+5,j+Npad+5);
        end
    end
    %draw the snapshot
  
if rem(n,isnap)==0
      imagesc(temp);
      axis equal
      axis([1 NX 1 NY]);
      colormap gray
      xlabel('x'),ylabel('z')
      title(sprintf('The current Timestep is %i ',n));
      drawnow
end
end
toc;

