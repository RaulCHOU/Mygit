clc
clear all
close all
N=1000;          %���������
d=0.012;      %���뻡��
delta1=0.3;   %�������1 
delta2=0.5;   %�������2
delta3=0.7;   %�������3
R0=35;        %��뾶
theta=d;      %��λ��
phi=asin(R0*sin(d/2)*sqrt(3)/R0); %����
A=[0,0];      %ԭ������
B=[theta/2,A(1,2)+phi];
C=[theta,A(1,2)];
Theta=[A(1,1),B(1,1),C(1,1)];
Phi=[A(1,2),B(1,2),C(1,2)];
eT=0;
eP=0;         %������
figure(1)
for n=1:N
    R=[rand()*theta,rand()*phi];   %���������ھ��ȷֲ�
%   R=[abs(0.2*randn(1,1)*theta),abs(0.2*randn(1,1)*phi)];   %���������ڷ��Ӿ�ֵΪ0,����Ϊ0.2db����̬�ֲ�
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
    plot(R(1,1),R(1,2),'*');  %�����
    hold on
    
%     Ee1=normrnd(0,sqrt(delta1/3));    %�÷������delta1����3����Χ������ֲ�Ҫ��������
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
    
    E=[1-(R(1,2)-A(1,2))/(2*phi)-(R(1,1)-A(1,1))/theta,(R(1,2)-A(1,2))/phi,(R(1,1)-A(1,1))/theta-(R(1,2)-A(1,2))/(2*phi)]; %���������߸���ֵ ��һ��
    E1=E;
%     E1=[E(1)+Ee1,E(2)+Ee2,E(3)+Ee3];   %������Ӿ�ֵΪ0������Ϊ1/3���������ķ������
%   E2=[E(1)+sqrt(delta2/3)*randn(1,1),E(2)+sqrt(delta2/3)*randn(1,1),E(3)+sqrt(delta2/3)*randn(1,1)];   %������Ӿ�ֵΪ0������Ϊ1/3���������ķ������
%   E3=[E(1)+sqrt(delta3/3)*randn(1,1),E(2)+sqrt(delta3/3)*randn(1,1),E(3)+sqrt(delta3/3)*randn(1,1)];   %������Ӿ�ֵΪ0������Ϊ1/3���������ķ������
    mole1=0;  %����1
    mole2=0;  %����2
    deno=0;   %��ĸ
    for i=1:3      
        mole1=E1(i)*Theta(i)+mole1;
        mole2=E1(i)*Phi(i)+mole2;
        deno=E1(i)+deno;
    end
    THETA=mole1/deno;   %�ϳ�Ŀ�귽λ��
    PHI=mole2/deno;     %�ϳ�Ŀ������
    
    e1(n)=R(1,1)-THETA;
    e2(n)=R(1,2)-PHI;
    
    % ��λ�����
    % ��ĸ
%     mm=E1(1,1)^2+E1(1,2)^2+E1(1,3)^2+2*(E1(1,1)*E1(1,2)+E1(1,2)*E1(1,3)+E1(1,3)*E1(1,1));
%     % ����
%     z=E1(1,1)^2*Theta(1,1)+E1(1,1)*(E1(1,2)*Theta(1,1)+E1(1,3)*Theta(1,1)+E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,3))+E1(1,2)^2*Theta(1,2)+E1(1,3)^2*Theta(1,3)+E1(1,2)*E1(1,3)*Theta(1,2)+E1(1,2)*E1(1,3)*Theta(1,3);
%     % ���Ӷ�E1��
%     zz1=(2*E1(1,1)*Theta(1,1)+E1(1,2)*Theta(1,1)+E1(1,3)*Theta(1,1)+E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,3))*mm-z*(2*E1(1,1)+2*E1(1,2)+2*E1(1,3));
%     % ����˸�������E1��ƫ��*E1���
%     EE1=zz1*Ee1/(mm^2);
%     zz2=(E1(1,1)*Theta(1,1)+E1(1,1)*Theta(1,2)+2*E1(1,2)*Theta(1,2)+E1(1,3)*Theta(1,2)+E1(1,3)*Theta(1,3))*mm-z*(2*E1(1,2)+2*E1(1,1)+2*E1(1,3));
%     EE2=zz2*Ee2/(mm^2);
%     zz3=(E1(1,1)*Theta(1,1)+E1(1,1)*Theta(1,3)+2*E1(1,3)*Theta(1,3)+E1(1,2)*Theta(1,2)+E1(1,2)*Theta(1,3))*mm-z*(2*E1(1,3)+2*E1(1,2)+2*E1(1,1));
%     EE3=zz3*Ee3/(mm^2);
%     ET(n)=EE1+EE2+EE3;  %��λ�����
%     
%     % ���������
%     z=E1(1,1)^2*Phi(1,1)+E1(1,1)*(E1(1,2)*Phi(1,1)+E1(1,3)*Phi(1,1)+E1(1,2)*Phi(1,2)+E1(1,3)*Phi(1,3))+E1(1,2)^2*Phi(1,2)+E1(1,3)^2*Phi(1,3)+E1(1,3)*Phi(1,2)+E1(1,2)*E1(1,3)*Phi(1,3);
%     % ���Ӷ�E1��
%     zz1=(2*E1(1,1)*Phi(1,1)+E1(1,2)*Phi(1,1)+E1(1,3)*Phi(1,1)+E1(1,2)*Phi(1,2)+E1(1,3)*Phi(1,3))*mm-z*(2*E1(1,1)+2*E1(1,2)+2*E1(1,3));
%     % ����˸�������E1��ƫ��*E1���
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
% countT=hist(ET,resultT);    %��������
% figure(2)
% plot(resultT,countT,'.');
% ylabel('��������');
% title('��������Ŀ��λ�÷�λ��λ�õ�Ӱ��');
% 
% resultP=unique(EP);
% countP=hist(EP,resultP);    %��������
% figure(3)
% plot(resultP,countP,'.');
% ylabel('��������');
% title('��������Ŀ��λ������λ�õ�Ӱ��');

%% �������״�Ͳ���ͼ
%% �ϳ�Ŀ���Դ
theta_a=0.01;   %��λ������ȣ����ȣ�
theta_b=0.01;   %����������ȣ����ȣ�
delta_a=R(1);       %�����ǲ� Ŀ��-����      
delta_b=R(2);       
phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
ka=1;
kb=1;               %��һ�������б��
m=[1,-1,-1,1];
n=[1,1,-1,-1];      %����ϵ��
xx=1.5708*(sin(1.18964*delta_a)-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
yy=1.5708*(sin(1.18964*delta_b)-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
F_sigma=0;      %�ͷ���ͼ
F_a=0;          %��λ���ͼ
F_b=0;          %�������ͼ 
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

%% ��Ԫ��Ϊ��Դ
theta_a=0.01;   %��λ������ȣ����ȣ�
theta_b=0.01;   %����������ȣ����ȣ�
delta_a=[Theta(1),Theta(2),Theta(3)];      
delta_b=[Phi(1),Phi(2),Phi(3)];        
phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
ka=1;
kb=1;               %��һ�������б��
m=[1,-1,-1,1];
n=[1,1,-1,-1];      %����ϵ��
for ii=1:3          %��Ԫ���Դ
    xx=1.5708*(sin(1.18964*delta_a(ii))-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
    yy=1.5708*(sin(1.18964*delta_b(ii))-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
    F_sigma=0;      %�ͷ���ͼ
    F_a=0;          %��λ���ͼ
    F_b=0;          %�������ͼ 
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

%% �ϳ�Ŀ���Դ
% k=0.73;
% d=0.16;
% labda=2*pi/k;
% theta_a=1.2*labda/d;   %��λ������ȣ����ȣ�
% theta_b=1.2*labda/d;   %����������ȣ����ȣ�
% delta_a=R(1);       %�����ǲ� Ŀ��-����      
% delta_b=R(2);       
% phi_i=[pi/4,3*pi/4,5*pi/4,7*pi/4];
% ka=0.707;
% kb=0.707    ;               %��һ�������б��
% m=[1,-1,-1,1];
% n=[1,1,-1,-1];      %����ϵ��
% xx=1.5708*(sin(1.18964*delta_a)-sin(theta_a*cos(phi_i)/sqrt(2)))/sin(theta_a/2);
% yy=1.5708*(sin(1.18964*delta_b)-sin(theta_b*cos(phi_i)/sqrt(2)))/sin(theta_b/2);
% F_sigma=0;      %�ͷ���ͼ
% F_a=0;          %��λ���ͼ
% F_b=0;          %�������ͼ 
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

%% ���ݺͲ���ͼ����E1��E2��E3
% x=0;    %����ͷ����������һ���ֱ������ϵ����
% y=0;
% z=0;
% Vs=1;   %��Ԫ��Ȩ����
% k=2*pi/1.5;    
% j=sqrt(-1);
% for i=1:3
%     Ri(i)=sqrt((x-R0*cos(Phi(i))*sin(Theta(i)))^2+(y-R0*cos(Phi(i))*cos(Theta(i)))^2+(z-R0*sin(Phi(i)))^2);
%     ti(i)=x*cos(Phi(i))*sin(Theta(i))+y*cos(Phi(i))*sin(Theta(i))+z*sin(Phi(i));
%     R=x^2+y^2+z^2;
%     Ri(i)=R0-ti(i)+1/2*(R-ti(i)^2)/R0;   %����ͷ��������һ��x y z����Ԫ�����Ԫ�ľ���
% 
%     Ei(i)=1;   %Դ����ֵ
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



