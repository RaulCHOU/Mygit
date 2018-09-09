clc
clear all
close all
N=1000;          %随机数个数
d=0.012;      %距离弧度
delta1=0.3;   %幅度误差1 
delta2=0.5;   %幅度误差2
delta3=0.7;   %幅度误差3
R0=35;        %球半径
theta=d;      %方位角
phi=asin(R0*sin(d/2)*sqrt(3)/R0); %仰角
A=[0,0];      %原点天线
B=[theta/2,A(1,2)+phi];
C=[theta,A(1,2)];
Theta=[A(1,1),B(1,1),C(1,1)];
Phi=[A(1,2),B(1,2),C(1,2)];
eT=0;
eP=0;         %求误差和
figure(1)
for n=1:N
    R=[rand()*theta,rand()*phi];   %在三角形内均匀分布
%   R=[abs(0.2*randn(1,1)*theta),abs(0.2*randn(1,1)*phi)];   %在三角形内服从均值为0,方差为0.2db的正态分布
    PA=[A(1,1)-R(1,1),A(1,2)-R(1,2)];
    PB=[C(1,1)-R(1,1),C(1,2)-R(1,2)];
    PC=[B(1,1)-R(1,1),B(1,2)-R(1,2)];
    t1=PA(1,1)*PC(1,2)-PA(1,2)*PC(1,1);
    t2=PC(1,1)*PB(1,2)-PC(1,2)*PB(1,1);
    t3=PB(1,1)*PA(1,2)-PB(1,2)*PA(1,1); 
    if ((t1*t2>=0 && t1*t3>=0)==0)
        if(R(1,1)>=B(1,1))
            R(1,1)=B(1,1)+C(1,1)-R(1,1);
            R(1,2)=B(1,2)+C(1,2)-R(1,2);
        else
            R(1,1)=A(1,1)+B(1,1)-R(1,1);
            R(1,2)=A(1,2)+B(1,2)-R(1,2);
        end
    end
    plot(R(1,1),R(1,2),'*');  %随机点
    hold on
    
%     Ee1=normrnd(0,sqrt(delta1/3));    %用幅度误差delta1产生3个范围内满足分布要求的随机数
%     while(abs(Ee1)>=delta1)
%         Ee1=normrnd(0,sqrt(delta1/3));
%     end
%     Ee2=normrnd(0,sqrt(delta1/3));
%     while(abs(Ee2) >= delta1)
%         Ee2=normrnd(0,sqrt(delta1/3));
%     end
%     Ee3=normrnd(0,sqrt(delta1/3));
%     while(abs(Ee3) >= delta1)
%         Ee3=normrnd(0,sqrt(delta1/3));
%     end
    
    E=[1-(R(1,2)-A(1,2))/(2*phi)-(R(1,1)-A(1,1))/theta,(R(1,2)-A(1,2))/phi,(R(1,1)-A(1,1))/theta-(R(1,2)-A(1,2))/(2*phi)]; %计算三天线各幅值 归一化
    E1=E;
%     E1=[E(1)+Ee1,E(2)+Ee2,E(3)+Ee3];   %引入服从均值为0，方差为1/3最大幅度误差的幅度误差
%   E2=[E(1)+sqrt(delta2/3)*randn(1,1),E(2)+sqrt(delta2/3)*randn(1,1),E(3)+sqrt(delta2/3)*randn(1,1)];   %引入服从均值为0，方差为1/3最大幅度误差的幅度误差
%   E3=[E(1)+sqrt(delta3/3)*randn(1,1),E(2)+sqrt(delta3/3)*randn(1,1),E(3)+sqrt(delta3/3)*randn(1,1)];   %引入服从均值为0，方差为1/3最大幅度误差的幅度误差
    mole1=0;  %分子1
    mole2=0;  %分子2
    deno=0;   %分母
    for i=1:3      
        mole1=E1(i)*Theta(i)+mole1;
        mole2=E1(i)*Phi(i)+mole2;
        deno=E1(i)+deno;
    end
    THETA=mole1/deno;   %合成目标方位角
    PHI=mole2/deno;     %合成目标仰角
    
    e1(n)=R(1,1)-THETA;
    e2(n)=R(1,2)-PHI;
    
    % 求方位角误差
    % 分母
%     mm=E1(1,1)^2+E1(1,2)^2+E1(1,3)^2+2*(E1(1,1)*E1(1,2)+E1(1,2)*E1(1,3)+E1(1,3)*E1(1,1));
%     % 分子
%     z=E1(1,1)^2*Theta(1,1)+E1(1,1)*(E1(1,2)*Theta(1,1)+E1(1,3)*Theta(1,1)+E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,3))+E1(1,2)^2*Theta(1,2)+E1(1,3)^2*Theta(1,3)+E1(1,2)*E1(1,3)*Theta(1,2)+E1(1,2)*E1(1,3)*Theta(1,3);
%     % 分子对E1求导
%     zz1=(2*E1(1,1)*Theta(1,1)+E1(1,2)*Theta(1,1)+E1(1,3)*Theta(1,1)+E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,3))*mm-z*(2*E1(1,1)+2*E1(1,2)+2*E1(1,3));
%     % 角闪烁方差对于E1的偏导*E1误差
%     EE1=zz1*Ee1/(mm^2);
%     zz2=(E1(1,1)*Theta(1,1)+E1(1,1)*Theta(1,2)+2*E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,2)+E1(1,3)*Theta(1,3))*mm-z*(2*E1(1,2)+2*E1(1,1)+2*E1(1,3));
%     EE2=zz2*Ee2/(mm^2);
%     zz3=(E1(1,1)*Theta(1,1)+E1(1,1)*Theta(1,3)+2*E1(1,3)*Theta(1,3)+E1(1,2)*Theta(1,2)+E1(1,2)*Theta(1,3))*mm-z*(2*E1(1,3)+2*E1(1,2)+2*E1(1,1));
%     EE3=zz3*Ee3/(mm^2);
%     ET(n)=EE1+EE2+EE3;  %方位角误差
%     
%     % 求仰角误差
%     z=E1(1,1)^2*Phi(1,1)+E1(1,1)*(E1(1,2)*Phi(1,1)+E1(1,3)*Phi(1,1)+E1(1,2)*Phi(1,2)+E1(1,3)*Phi(1,3))+E1(1,2)^2*Phi(1,2)+E1(1,3)^2*Phi(1,3)+E1(1,3)*Phi(1,2)+E1(1,2)*E1(1,3)*Phi(1,3);
%     % 分子对E1求导
%     zz1=(2*E1(1,1)*Phi(1,1)+E1(1,2)*Phi(1,1)+E1(1,3)*Phi(1,1)+E1(1,2)*Phi(1,2)+E1(1,3)*Phi(1,3))*mm-z*(2*E1(1,1)+2*E1(1,2)+2*E1(1,3));
%     % 角闪烁方差对于E1的偏导*E1误差
%     EE1=zz1*Ee1/(mm^2);
%     zz2=(E1(1,1)*Phi(1,1)+E1(1,1)*Phi(1,2)+2*E1(1,2)*Phi(1,2)+E1(1,3)*Phi(1,2)+E1(1,3)*Phi(1,3))*mm-z*(2*E1(1,2)+2*E1(1,1)+2*E1(1,3));
%     EE2=zz2*Ee2/(mm^2);
%     zz3=(E1(1,1)*Phi(1,1)+E1(1,1)*Phi(1,3)+2*E1(1,3)*Phi(1,3)+E1(1,2)*Phi(1,2)+E1(1,2)*Phi(1,3))*mm-z*(2*E1(1,3)+2*E1(1,2)+2*E1(1,1));
%     EE3=zz3*Ee3/(mm^2);
%     EP(n)=EE1+EE2+EE3;
%     
%     eT=ET(n)+eT;
%     eP=EP(n)+eP;
end
% eT=eT/n
% eP=eP/n
% for i=1:N
%     varT=(ET(i)-eT)^2;
%     varP=(EP(i)-eP)^2;
% end
% varT=varT/i
% varP=varP/i

% resultT=unique(ET);
% countT=hist(ET,resultT);    %按误差计数
% figure(2)
% plot(resultT,countT,'.');
% ylabel('发生次数');
% title('幅度误差对目标位置方位角位置的影响');
% 
% resultP=unique(EP);
% countP=hist(EP,resultP);    %按误差计数
% figure(3)
% plot(resultP,countP,'.');
% ylabel('发生次数');
% title('幅度误差对目标位置仰角位置的影响');

%% 求单脉冲雷达和差差方向图
%% 合成目标点源
theta_a=0.01;   %方位波束宽度（弧度）
theta_b=0.01;   %俯仰波束宽度（弧度）
delta_a=R(1);       %俯仰角差 目标-天线      
delta_b=R(2);       
phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
ka=1;
kb=1;               %归一化角误差斜率
m=[1,-1,-1,1];
n=[1,1,-1,-1];      %象限系数
xx=1.5708*(sin(1.18964*delta_a)-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
yy=1.5708*(sin(1.18964*delta_b)-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
F_sigma=0;      %和方向图
F_a=0;          %方位差方向图
F_b=0;          %俯仰差方向图 
for i=1:4
    f=sin(xx(i))*sin(yy(i))/(xx(i)*yy(i));
    f_sigma(i)=0.616853*f;
    f_a(i)=0.2595*ka*m(i)*f;
    f_b(i)=0.2595*kb*n(i)*f;
    F_sigma=f_sigma(i)+F_sigma;
    F_a=f_a(i)+F_a;
    F_b=f_b(i)+F_b;
end
Sum=F_sigma;
Daz=F_a;
Del=F_b;

%% 三元组为点源
theta_a=0.01;   %方位波束宽度（弧度）
theta_b=0.01;   %俯仰波束宽度（弧度）
delta_a=[Theta(1),Theta(2),Theta(3)];      
delta_b=[Phi(1),Phi(2),Phi(3)];        
phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
ka=1;
kb=1;               %归一化角误差斜率
m=[1,-1,-1,1];
n=[1,1,-1,-1];      %象限系数
for ii=1:3          %三元组点源
    xx=1.5708*(sin(1.18964*delta_a(ii))-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
    yy=1.5708*(sin(1.18964*delta_b(ii))-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
    F_sigma=0;      %和方向图
    F_a=0;          %方位差方向图
    F_b=0;          %俯仰差方向图 
    for i=1:4
        f=sin(xx(i))*sin(yy(i))/(xx(i)*yy(i));
        f_sigma(i)=0.616853*f;
        f_a(i)=0.2595*ka*m(i)*f;
        f_b(i)=0.2595*kb*n(i)*f;
        F_sigma=f_sigma(i)+F_sigma;
        F_a=f_a(i)+F_a;
        F_b=f_b(i)+F_b;
    end
    SUM(ii)=F_sigma;
    DAZ(ii)=F_a;
    DEL(ii)=F_b;
end
syms E1 E2 E3;
[E1,E2,E3]=solve('DAZ(1)*E1+DAZ(2)*E2+DAZ(3)*E3=0',...
    'DEL(1)*E1+DEL(2)*E2+DEL(3)*E3=0',...
    'E1+E2+E3=1');
E1=eval(E1)
E2=eval(E2)
E3=eval(E3)

%% 合成目标点源
% k=0.73;
% d=0.16;
% labda=2*pi/k;
% theta_a=1.2*labda/d;   %方位波束宽度（弧度）
% theta_b=1.2*labda/d;   %俯仰波束宽度（弧度）
% delta_a=R(1);       %俯仰角差 目标-天线      
% delta_b=R(2);       
% phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
% ka=0.707;
% kb=0.707    ;               %归一化角误差斜率
% m=[1,-1,-1,1];
% n=[1,1,-1,-1];      %象限系数
% xx=1.5708*(sin(1.18964*delta_a)-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
% yy=1.5708*(sin(1.18964*delta_b)-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
% F_sigma=0;      %和方向图
% F_a=0;          %方位差方向图
% F_b=0;          %俯仰差方向图 
% for i=1:4
%     f=sin(xx(i))*sin(yy(i))/(xx(i)*yy(i));
%     f_sigma(i)=0.616853*f;
%     f_a(i)=0.2595*ka*m(i)*f;
%     f_b(i)=0.2595*kb*n(i)*f;
%     F_sigma=f_sigma(i)+F_sigma;
%     F_a=f_a(i)+F_a;
%     F_b=f_b(i)+F_b;
% end
% SUM=F_sigma;
% DAZ=F_a;
% DEL=F_b;
% syms E1 E2 E3;
% [E1,E2,E3]=solve('abs(DAZ(1)*abs(E1)+DAZ(2)*abs(E2)+DAZ(3)*abs(E3))=0',...
%     'abs(DEL(1)*abs(E1)+DEL(2)*abs(E2)+DEL(3)*abs(E3))=0',...
%     'abs(abs(E1)+abs(E2)+abs(E3)-1)=0');
% E1
% E2
% E3

%% 根据和差差方向图反馈E1、E2、E3
% x=0;    %导引头天线上任意一点的直角坐标系坐标
% y=0;
% z=0;
% Vs=1;   %阵元加权因子
% k=2*pi/1.5;    
% j=sqrt(-1);
% for i=1:3
%     Ri(i)=sqrt((x-R0*cos(Phi(i))*sin(Theta(i)))^2+(y-R0*cos(Phi(i))*cos(Theta(i)))^2+(z-R0*sin(Phi(i)))^2);
%     ti(i)=x*cos(Phi(i))*sin(Theta(i))+y*cos(Phi(i))*sin(Theta(i))+z*sin(Phi(i));
%     R=x^2+y^2+z^2;
%     Ri(i)=R0-ti(i)+1/2*(R-ti(i)^2)/R0;   %导引头天线上任一点x y z到三元组各阵元的距离
% 
%     Ei(i)=1;   %源幅度值
%     fo(i)=sqrt(1-cos(Phi(i))^2*sin(R(1,1)-Theta(1))^2);
%     ph(i)=Vs*fo(i)*exp(j*k*(ti(i)*(R-ti(i)^2)/(-2*R0)));
%    
% end
% F=E1(1)*ph(1)+E1(2)*ph(2)+E1(3)*ph(3);
% F1=4*F;
% F2=2*F;
% F3=2*F;
% F4=4*F;
% SUM=F1+F2+F3+F4;
% DDAZ=(F1+F4-F2-F3)/SUM;
% DDEL=(F1+F2-F3-F4)/SUM;
% for i=1:3
%     DAZ(i)=4*ph(i)/SUM;
%     DEL(i)=0;
% end
% syms E1 E2 E3;
% [E1,E2,E3]=solve('abs(DAZ(1)*abs(E1)+DAZ(2)*abs(E2)+DAZ(3)*abs(E3))=0',...
%   'abs(DEL(1)*abs(E1)+DEL(2)*abs(E2)+DEL(3)*abs(E3))=0',...
%   'abs(abs(E1)+abs(E2)+abs(E3)-1)=0');
%     



