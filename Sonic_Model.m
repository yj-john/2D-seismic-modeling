%2D acoustic modeling
%规则网格正演模拟
%2阶时间差分精度，10阶空间差分精度
clc;
clear;
tic;
X=3000,Y=3000;%模型计算区域大小
DX=10,DY=10;%模型的网格间距
nx=X/DX+10,ny=Y/DY+10;%模型网格扩充大小
NX=X/DX,NY=Y/DY;%模型实际大小
SX=1500,SY=1500;%震源位置
sx=SX/DX+5,sy=SY/DY+5;%震源网格位置
Vp=1500*ones(nx,ny);
nt=0.8;%正演模拟时间
dt=0.001;%时间步长间隔
ft=40;%震源主频大小
%雷克子波
ns=100%子波的采样点数
it=(1:100)*dt;
wave=(1-2*(pi*ft.*(it-1/ft).^2)).*exp(-(pi*ft.*(it-1/ft)).^2);
figure(1); 
plot(1:ns,wave);xlabel('采样间隔');
%声波正演运行
%差分系数大小
fdcoeff=[-2.9672,1.6666667,-0.2380952,0.0396825,-0.0049603,0.0003175];
isnap=5;
%波场初始化
p=zeros(nx,ny);%当前时候的波场
p1=zeros(nx,ny);%前一个时刻的波场
p2=zeros(nx,ny);%后一个时刻的波场
d2px=zeros(nx,ny);
d2py=zeros(nx,ny);
%波场快照初始化
temp=zeros(NX,NY);
tdsnapshot=zeros(NX,NY,nt/dt);
%波场延拓
figure(2);
for n=1:nt/dt
  for i=6:nx-5
     for j=6:ny-5
           d2px(i,j)=fdcoeff(1)*p(i,j);
           d2py(i,j)=fdcoeff(1)*p(i,j);
        for m=2:length(fdcoeff)
           d2px(i,j)=d2px(i,j)+fdcoeff(m)*(p(i+m-1,j)+p(i-m+1,j));
           d2py(i,j)=d2py(i,j)+fdcoeff(m)*(p(i,j+m-1)+p(i,j-m+1));
        end
         p2(i,j)=2*p(i,j)-p1(i,j)+Vp(i,j)^2*(dt^2)*(d2px(i,j)/DX^2+d2py(i,j)/DY^2);
         
     end
  end
%震源函数加载
if n<ns
    p2(sx,sy)=wave(n);
end
%波场更新
p1=p;
p=p2;
%-----------------------------计算波场快照------------------------------%
% 输出波场快照
for i=1:NX
    for j=1:NY
     temp(i,j)=p(i+5,j+5);
     tdsnapshot(i,j,n)=p(i+5,j+5);
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
      
