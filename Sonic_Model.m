%2D acoustic modeling
%������������ģ��
%2��ʱ���־��ȣ�10�׿ռ��־���
clc;
clear;
tic;
X=3000,Y=3000;%ģ�ͼ��������С
DX=10,DY=10;%ģ�͵�������
nx=X/DX+10,ny=Y/DY+10;%ģ�����������С
NX=X/DX,NY=Y/DY;%ģ��ʵ�ʴ�С
SX=1500,SY=1500;%��Դλ��
sx=SX/DX+5,sy=SY/DY+5;%��Դ����λ��
Vp=1500*ones(nx,ny);
nt=0.8;%����ģ��ʱ��
dt=0.001;%ʱ�䲽�����
ft=40;%��Դ��Ƶ��С
%�׿��Ӳ�
ns=100%�Ӳ��Ĳ�������
it=(1:100)*dt;
wave=(1-2*(pi*ft.*(it-1/ft).^2)).*exp(-(pi*ft.*(it-1/ft)).^2);
figure(1); 
plot(1:ns,wave);xlabel('�������');
%������������
%���ϵ����С
fdcoeff=[-2.9672,1.6666667,-0.2380952,0.0396825,-0.0049603,0.0003175];
isnap=5;
%������ʼ��
p=zeros(nx,ny);%��ǰʱ��Ĳ���
p1=zeros(nx,ny);%ǰһ��ʱ�̵Ĳ���
p2=zeros(nx,ny);%��һ��ʱ�̵Ĳ���
d2px=zeros(nx,ny);
d2py=zeros(nx,ny);
%�������ճ�ʼ��
temp=zeros(NX,NY);
tdsnapshot=zeros(NX,NY,nt/dt);
%��������
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
%��Դ��������
if n<ns
    p2(sx,sy)=wave(n);
end
%��������
p1=p;
p=p2;
%-----------------------------���㲨������------------------------------%
% �����������
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
      
