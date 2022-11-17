%2D acoustic modeling using staggered grid
%交错网格声波正演模拟
%十阶交错网格精度
clc;
clear;
tic;
X=1000,Y=1000;%模型计算区域大小
DX=10,DY=10;%网格间距大小
Npad=20;
NX=X/DX,NY=Y/DY;%网格实际数
nx=NX+10+2*Npad,ny=NY+10+2*Npad;%网格扩充后的大小
SX=500,SY=500;%震源位置
sx=SX/DX+5+Npad,sy=SY/DY+5+Npad;%%震源网格位置
Vp=1500*ones(nx,ny);
e=2*ones(nx,ny);
nt=0.6%正演模拟时间
dt=0.001%采样时间间隔
ft=40%主频大小
%雷克子波
ns=100;
it=(1:100)*dt;
wave=(1-2*(pi*ft.*(it-1/ft).^2)).*exp(-(pi*ft.*(it-1/ft)).^2);
figure(1); 
plot(1:ns,wave);xlabel('采样间隔');
%声波正演运行
%交错网格差分系数
fdcoeff=[1.1211  -0.089722 0.013843 -0.0017657 0.00011868]
isnap=1;
%波场初始化
vx1=zeros(nx,ny);%前一个时刻(k)的Vx分量
vx2=zeros(nx,ny);%后一个时刻(k-1)的Vx分量
vy1=zeros(nx,ny);%前一个时刻(k)的Vy分量
vy2=zeros(nx,ny);%后一个时刻(k-1)的Vy分量
p=zeros(nx,ny);
p1=zeros(nx,ny);
px=zeros(nx,ny);%前一个时刻(k+1/2)的分量
py=zeros(nx,ny);
px1=zeros(nx,ny);%后一个时刻(k-1/2)的分量
py1=zeros(nx,ny);
dx=zeros(nx,ny);
dy=zeros(nx,ny);
%波场快照初始化
temp=zeros(NX,NY);
%吸收边界系数初始化
 %完全匹配层  其衰减模型公式d1(x)=3*Vmax/(2*npad)*log(1/R)(x/npad)^2
 %d2(z)=3*Vmax/(2*npad)*log(1/R)(z/npad)^2
 % left //左边MPL边界衰减
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
% up 向上边界衰减
   for i=Npad+6:(nx-(Npad+5))
       for j=6:Npad+5
           dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(Npad+6-j)/Npad]^2;
           dx(i,j)=0;
       end
   end
%down 向下边界衰减
   for i=Npad+6:(nx-(Npad+5))
      for j=(ny-Npad-4):(ny-5)
            dy(i,j)=3*Vp(i,j)/(2*Npad)*log(10^6)*[(j-(ny-5-Npad))/Npad]^2;
            dx(i,j)=0;
      end
   end
%四个边角
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
%波场延拓
for n=1:nt/dt
dvx=zeros(nx,ny);%Vx分量导数
dvy=zeros(nx,ny);%Vy分量导数
dpx=zeros(nx,ny);%dp的x方向导数
dpy=zeros(nx,ny);%dp的y方向导数

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
 %震源函数加载
    if n<ns
        p(sx,sy)=wave(n);
    end
      %波场更新
    p1=p;
    px1=px;
    py1=py;
   vx2=vx1;
   vy2=vy1;
    %输出波场快照
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

