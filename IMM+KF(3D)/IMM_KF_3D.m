clear all
clc
echo off
tic

%==========================================================================
%IMM+KF 3D
%==========================================================================

%�������
simTime = 100;%�����������
T = 0.05;%����ʱ��
G = 1;%����������������Ϊ�˷���ֱ����1��
H = [1,0,0,0,0,0,0,0,0;
     0,0,0,1,0,0,0,0,0;
     0,0,0,0,0,0,1,0,0];
r = 20;%������������
q = 5;%������������
R = r*eye(3);%ģ����������Э�������
Q = q*eye(9);%ģ�͹�������Э�������
R1 = sqrtm(R)*randn(3,simTime);
Q1 = sqrtm(Q)*randn(9,simTime);
w1 = 3*2*pi/360;
w2 = -3*2*pi/360;
w3 = 6*2*pi/360;%ת����-3��

%����ģ��(cv)
F_CV = zeros(9,9);
F1 = [1 T 0;
      0 1 0;
      0 0 0];
F_CV(1:3,1:3) = F1;F_CV(4:6,4:6) = F1;F_CV(7:9,7:9) = F1;
%�ȼ���ģ��(ca)
F_CA = zeros(9,9);
F2 = [1 T T^2/2;
      0 1 T;
      0 1 1];
F_CA(1:3,1:3) = F2;F_CA(4:6,4:6) = F2;F_CA(7:9,7:9) = F2;
%������Эͬת��ģ��(csct)
F_CSCT = zeros(9,9);
F3_1 = [1 sin(w1*T)/w1 (1-cos(w1*T))/w1^2;
      0 cos(w1*T) sin(w1*T)/w1;
      0 -w1*sin(w1*T) cos(w1*T)];
F3_2 = [1 sin(w2*T)/w1 (1-cos(w2*T))/w2^2;
      0 cos(w2*T) sin(w2*T)/w2;
      0 -w2*sin(w2*T) cos(w2*T)];
F3_3 = [1 sin(w3*T)/w3 (1-cos(w3*T))/w3^2;
      0 cos(w3*T) sin(w3*T)/w3;
      0 -w3*sin(w3*T) cos(w3*T)];
F_CSCT(1:3,1:3) = F3_1;F_CSCT(4:6,4:6) = F3_2;F_CSCT(7:9,7:9) = F3_3;

%������������
x0 = [20,20,0,20,20,0,20,20,0]';
rand('state',sum(100*clock));
x = zeros(9,simTime);%����Ŀ��״̬
x_true = zeros(9,simTime);%��ʵĿ��״̬
z = zeros(3,simTime);%��ʵĿ�����ֵ
z_true = zeros(3,simTime);%ʵ��Ŀ�����ֵ
x(:,1) = x0;
x_true(:,1) = x0;
z(:,1) = H*x(:,1)+R1(:,1);%�������������������ݣ���ʼ��
z_true(:,1) = H*x(:,1);%������������ʵ�������ݣ���ʼ��

%��ʵ����ֵ��Ԥ��ֵ
for a = 2:simTime
    if (a>=20)&&(a<=40)
         x_true(:,a) = F_CA*x_true(:,a-1);
    elseif (a>=60)&&(a<=80)
         x_true(:,a) = F_CSCT*x_true(:,a-1);
    else
         x_true(:,a) = F_CV*x_true(:,a-1);
    end
    
    z_true(:,a) = H*x_true(:,a);
    
end
%ʵ�ʲ���ֵ��Ԥ��ֵ
for a = 2:simTime
    if (a>=20)&&(a<=40)
         x(:,a) = F_CA*x_true(:,a-1)+G*Q1(:,a-1);
    elseif (a>=60)&&(a<=80)
         x(:,a) = F_CSCT*x_true(:,a-1)+G*Q1(:,a-1);
    else
         x(:,a) = F_CV*x_true(:,a-1)+G*Q1(:,a-1);
    end
    
    z(:,a) = H*x(:,a)+R1(:,a);
    
end

%==========================================================================
%IMM
%==========================================================================
%��һ�� ��ʼ��
x1_IMM = zeros(9,1);%ģ��1 IMM�㷨״̬����ֵ
x2_IMM = zeros(9,1);
x3_IMM = zeros(9,1);
x_pro_IMM = zeros(9,simTime);%IMM�㷨�ۺ�״̬����ֵ
P_IMM = zeros(9,9,simTime);%IMM�㷨ģ���ۺ�״̬Э�������
P1_IMM = zeros(9,9);%ģ��1״̬Э�������
P2_IMM = zeros(9,9);
P3_IMM = zeros(9,9);
r1_IMM = zeros(3,1);%ģ��1Ԥ��Ĺ۲�ֵ��۲�ֵ�Ĳв�
r2_IMM = zeros(3,1);
r3_IMM = zeros(3,1);
S1_IMM = zeros(3,3);%ģ��1Ԥ��Ĺ۲�ֵ��۲�ֵ�Ĳв��Э�������
S2_IMM = zeros(3,3);
S3_IMM = zeros(3,3);

x_pro_IMM(:,1) = x0;
pij = [0.8,0.1,0.1;
       0.1,0.8,0.1;
       0.1,0.1,0.8];%ģ�͸���ת�ƾ���
u_IMM = zeros(3,simTime);%IMM�㷨ģ�͸���
u_IMM(:,1) = [0.3,0.3,0.4];%IMM�㷨ģ�͸��ʳ�ʼֵ

x1_IMM = x0; x2_IMM = x0; x3_IMM = x0;

P0 = diag([0,0,0,0,0,0,0,0,0]);%��ʼЭ�������
P1_IMM = P0; P2_IMM = P0; P3_IMM = P0;

for t = 1:simTime-1
    %t = 1;
    
    %��һ�� ���뽻��
    c_j = pij'*u_IMM(:,t);
    
    ui1 = (1/c_j(1))*pij(:,1).*u_IMM(:,t);%����ģ�ͻ�ϸ���
    ui2 = (1/c_j(2))*pij(:,2).*u_IMM(:,t);
    ui3 = (1/c_j(3))*pij(:,3).*u_IMM(:,t);
    
    %�����ģ�ͳ�ʼ���˲���״̬��������
    x01 = x1_IMM*ui1(1)+x2_IMM*ui1(2)+x3_IMM*ui1(3);%�����ģ�ͳ�ʼ״̬
    x02 = x1_IMM*ui2(1)+x2_IMM*ui2(2)+x3_IMM*ui2(3);
    x03 = x1_IMM*ui3(1)+x2_IMM*ui3(2)+x3_IMM*ui3(3);
    
    P01 = (P1_IMM+(x1_IMM-x01)*(x1_IMM-x01)')*ui1(1)+...%�����ģ�ͳ�ʼ״̬Э�������
          (P2_IMM+(x2_IMM-x01)*(x2_IMM-x01)')*ui1(2)+...
          (P3_IMM+(x3_IMM-x01)*(x3_IMM-x01)')*ui1(3);
    P02 = (P1_IMM+(x1_IMM-x02)*(x1_IMM-x02)')*ui2(1)+...
          (P2_IMM+(x2_IMM-x02)*(x2_IMM-x02)')*ui2(2)+...
          (P3_IMM+(x3_IMM-x02)*(x3_IMM-x02)')*ui2(3);
    P03 = (P1_IMM+(x1_IMM-x03)*(x1_IMM-x03)')*ui3(1)+...
          (P2_IMM+(x2_IMM-x03)*(x2_IMM-x03)')*ui3(2)+...
          (P3_IMM+(x3_IMM-x03)*(x3_IMM-x03)')*ui3(3);
      
    %�ڶ��� �������˲� 
    [x1_IMM,P1_IMM,r1_IMM,S1_IMM] = Kalman(x01,P01,z(:,t+1),F_CV,G,Q,H,R);
    [x2_IMM,P2_IMM,r2_IMM,S2_IMM] = Kalman(x02,P02,z(:,t+1),F_CA,G,Q,H,R);
    [x3_IMM,P3_IMM,r3_IMM,S3_IMM] = Kalman(x03,P03,z(:,t+1),F_CSCT,G,Q,H,R);

%     %�ڶ��� �޼��������˲�
%     [x1_IMM,P1_IMM,r1_IMM,S1_IMM] = U_Kalman(x01,P01,z(:,t+1),F1,G,Q,H,R);
%     [x2_IMM,P2_IMM,r2_IMM,S2_IMM] = U_Kalman(x02,P02,z(:,t+1),F2,G,Q,H,R);
%     [x3_IMM,P3_IMM,r3_IMM,S3_IMM] = U_Kalman(x03,P03,z(:,t+1),F3,G,Q,H,R); 
    
    %������ ģ�͸��ʸ���
    [u_IMM(:,t+1)] = Model_P_up(r1_IMM,r2_IMM,r3_IMM,S1_IMM,S2_IMM,S3_IMM,c_j);
    
    %���Ĳ� ģ������ۺ�
    [x_pro_IMM(:,t+1),P_IMM(:,:,t+1)] = Model_mix(x1_IMM,x2_IMM,x3_IMM,P1_IMM,P2_IMM,P3_IMM,u_IMM(:,t+1));
    
    
end
toc

%==========================================================================
%������������
%==========================================================================
z_IMMKF = H*x_pro_IMM(:,:); 
Rmse_m = Compute_Rmse_z(z,z_true,simTime);%�������ֵ�ľ��������
Rmse_IMMKF = Compute_Rmse_z(z_IMMKF,z_true,simTime);%�������ֵ�ľ��������

%==========================================================================
%��ͼ
%==========================================================================
%ģ�͸���
figure(1)
t  = 1:simTime;
plot(t,u_IMM(1,t),'k:',t,u_IMM(2,t),'r-',t,u_IMM(3,t),'b--','linewidth',2);grid on
title('IMM�㷨ģ�͸�������');
xlabel('t/s');ylabel('ģ�͸���');
legend('ģ��1','ģ��2','ģ��3');

figure(2)
plot3(z_true(1,:),z_true(2,:),z_true(3,:),'black');grid on;hold on;
plot3(z(1,:),z(2,:),z(3,:),'r');
plot3(z_IMMKF(1,:),z_IMMKF(2,:),z_IMMKF(3,:),'b');
hold off


magnify