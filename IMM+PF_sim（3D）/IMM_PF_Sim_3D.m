clear all
clc
echo off
tic

%==========================================================================
%IMM+PF 3D（simple）
%==========================================================================

%仿真参数
simTime = 100;%仿真迭代次数
T = 0.05;%采样时间
N = 50;%粒子滤波粒子数100
G = 1;%过程噪声驱动矩阵为了方便直接填1了
H = [1,0,0,0,0,0,0,0,0;
     0,0,0,1,0,0,0,0,0;
     0,0,0,0,0,0,1,0,0];
r = 20;%测量噪声方差
q = 5;%过程噪声方差
R = r*eye(3);%模型量测噪声协方差矩阵
Q = q*eye(9);%模型过程噪声协方差矩阵
R1 = sqrtm(R)*randn(3,simTime);
Q1 = sqrtm(Q)*randn(9,simTime);
w1 = 3*2*pi/360;
w2 = -3*2*pi/360;
w3 = 6*2*pi/360;%转弯率-3度

%匀速模型(cv)
F_CV = zeros(9,9);
F1 = [1 T 0;
      0 1 0;
      0 0 0];
F_CV(1:3,1:3) = F1;F_CV(4:6,4:6) = F1;F_CV(7:9,7:9) = F1;
%匀加速模型(ca)
F_CA = zeros(9,9);
F2 = [1 T T^2/2;
      0 1 T;
      0 1 1];
F_CA(1:3,1:3) = F2;F_CA(4:6,4:6) = F2;F_CA(7:9,7:9) = F2;
%常速率协同转弯模型(csct)
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

%产生测量数据
x0 = [20,20,0,20,20,0,20,20,0]';
rand('state',sum(100*clock));
x = zeros(9,simTime);%测量目标状态
x_true = zeros(9,simTime);%真实目标状态
z = zeros(3,simTime);%真实目标测量值
z_true = zeros(3,simTime);%实际目标测量值
x(:,1) = x0;
x_true(:,1) = x0;
z(:,1) = H*x(:,1)+R1(:,1);%产生含噪声的量测数据（初始）
z_true(:,1) = H*x(:,1);%不带噪声的真实量测数据（初始）

%真实测量值与预测值
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
%实际测量值与预测值
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
%第一步初始化
Particle_X1 = zeros(9,N);%从上一时刻采样的各模型的N个粒子状态
Particle_X2 = zeros(9,N);
Particle_X3 = zeros(9,N);
Particle_P1 = zeros(9,9,simTime);%这边协方差矩阵就不考虑所有粒子的了
Particle_P2 = zeros(9,9,simTime);
Particle_P3 = zeros(9,9,simTime);

x01 = zeros(9,N);%各模型交互后的状态
x02 = zeros(9,N);
x03 = zeros(9,N);
P01 = zeros(9,9,simTime);
P02 = zeros(9,9,simTime);
P03 = zeros(9,9,simTime);

x01_pf = zeros(9,N);%粒子滤波后各模型的N个粒子状态
x02_pf = zeros(9,N);
x03_pf = zeros(9,N);
P01_pf = zeros(size(P01));
P02_pf = zeros(size(P02));
P03_pf = zeros(size(P03));

x1_IMM = zeros(9,simTime);%各模型的输出期望
x2_IMM = zeros(9,simTime);
x3_IMM = zeros(9,simTime);
P1_IMM = zeros(9,9,simTime);%各模型状态误差的期望
P2_IMM = zeros(9,9,simTime);
P3_IMM = zeros(9,9,simTime);

r1_IMM = zeros(3,1);%这边也不考虑左右粒子的残差了
r2_IMM = zeros(3,1);
r3_IMM = zeros(3,1);

S1 = zeros(3,3);%各模型各粒子的残差的协方差矩阵
S2 = zeros(3,3);
S3 = zeros(3,3);

Lamda_1 = zeros(1,N);%模型1 似然函数
Lamda_2 = zeros(1,N);
Lamda_3 = zeros(1,N);
Lamda_11 = zeros(1,N);%模型1 似然函数
Lamda_21 = zeros(1,N);
Lamda_31 = zeros(1,N);

x_pro_IMMPF = zeros(9,simTime);%模型输出综合
x_pro_IMMPF_N = zeros(9,N);

P_IMMPF = zeros(9,9,simTime);%模型协方差矩阵综合

u_IMM_N = zeros(3,simTime);%模型更新概率

%初始化
pij = [0.8,0.1,0.1;
       0.1,0.8,0.1;
       0.1,0.1,0.8];%模型概率转移矩阵
   
u_IMM_N(:,1) = [0.3 0.3 0.4]';%设置初始模型概率

P0 = diag([10,10,10,10,10,10,10,10,10]);%设置初始协方差矩阵
x1_IMM(:,1) = x0;%设置初始各模型上一时刻输出期望
x2_IMM(:,1) = x0;
x3_IMM(:,1) = x0;
P1_IMM(:,:,1) = P0;
P2_IMM(:,:,1) = P0;
P3_IMM(:,:,1) = P0;
x_pro_IMMPF(:,1) = x0;


for t = 1:simTime-1
    %t = 1;
    
    %第一步 输入交互
    c_j = pij'*u_IMM_N(:,t);
    
    ui1 = (1./c_j(1))*pij(:,1).*u_IMM_N(:,t);
    ui2 = (1./c_j(2))*pij(:,2).*u_IMM_N(:,t);
    ui3 = (1./c_j(3))*pij(:,3).*u_IMM_N(:,t);
    
    %计算各模型初始化滤波的状态（交互）
    x01 = x1_IMM(:,t)*ui1(1)+x2_IMM(:,t)*ui1(2)+x3_IMM(:,t)*ui1(3);%计算各模型初始状态
    x02 = x1_IMM(:,t)*ui2(1)+x2_IMM(:,t)*ui2(2)+x3_IMM(:,t)*ui2(3);
    x03 = x1_IMM(:,t)*ui3(1)+x2_IMM(:,t)*ui3(2)+x3_IMM(:,t)*ui3(3);
    
    P01 = (P1_IMM(:,:,t)+(x1_IMM(:,t)-x01)*(x1_IMM(:,t)-x01)')*ui1(1)+...%计算各模型初始状态协方差矩阵
          (P2_IMM(:,:,t)+(x2_IMM(:,t)-x01)*(x2_IMM(:,t)-x01)')*ui1(2)+...
          (P3_IMM(:,:,t)+(x3_IMM(:,t)-x01)*(x3_IMM(:,t)-x01)')*ui1(3);
    P02 = (P1_IMM(:,:,t)+(x1_IMM(:,t)-x02)*(x1_IMM(:,t)-x02)')*ui2(1)+...
          (P2_IMM(:,:,t)+(x2_IMM(:,t)-x02)*(x2_IMM(:,t)-x02)')*ui2(2)+...
          (P3_IMM(:,:,t)+(x3_IMM(:,t)-x02)*(x3_IMM(:,t)-x02)')*ui2(3);
    P03 = (P1_IMM(:,:,t)+(x1_IMM(:,t)-x03)*(x1_IMM(:,t)-x03)')*ui3(1)+...
          (P2_IMM(:,:,t)+(x2_IMM(:,t)-x03)*(x2_IMM(:,t)-x03)')*ui3(2)+...
          (P3_IMM(:,:,t)+(x3_IMM(:,t)-x03)*(x3_IMM(:,t)-x03)')*ui3(3);
      
    %第二步 粒子滤波
    
    %从交互结果中抽取N个粒子 抽取范围为交互后的各模型协方差矩阵
    randn('state',sum(100*clock));
    QQ1 = sqrtm(P1_IMM(:,:,t))*randn(9,N);
    QQ2 = sqrtm(P2_IMM(:,:,t))*randn(9,N);
    QQ3 = sqrtm(P3_IMM(:,:,t))*randn(9,N);
    
    for i = 1:N
        
        Particle_X1(:,i) = x01+QQ1(:,i);%从上一时刻各模型粒子滤波的期望状态进行采样
        Particle_X2(:,i) = x02+QQ2(:,i);
        Particle_X3(:,i) = x03+QQ3(:,i);
        
    end
    
    %聚类
    [Index_cluster_X1,Particle_X1] = cluster(Particle_X1,x(:,t),5,N);
    [Index_cluster_X2,Particle_X2] = cluster(Particle_X2,x(:,t),5,N);
    [Index_cluster_X3,Particle_X3] = cluster(Particle_X3,x(:,t),5,N);
    
     %对每个模型分别进行滤波 得到 各粒子的权重 个粒子状态 与 各粒子与测量值之差ero（3维）
    [weight_model_1,x01_pf,ero1] = Particle_Residual_sim(N,Particle_X1,z(:,t+1),F_CV,H,r,q,0.0001);
    [weight_model_2,x02_pf,ero2] = Particle_Residual_sim(N,Particle_X2,z(:,t+1),F_CA,H,r,q,0.0001);
    [weight_model_3,x03_pf,ero3] = Particle_Residual_sim(N,Particle_X3,z(:,t+1),F_CSCT,H,r,q,0.0001);
    
    %对各模型的粒子求状态的均值 与 协方差矩阵的均值 与 残差协方差矩阵的均值
    x1_IMM(:,t+1) = [mean(x01_pf(1,:)) mean(x01_pf(2,:)) mean(x01_pf(3,:)) mean(x01_pf(4,:)) ...
                     mean(x01_pf(5,:)) mean(x01_pf(6,:)) mean(x01_pf(7,:)) mean(x01_pf(8,:)) ...
                     mean(x01_pf(9,:))]';
    x2_IMM(:,t+1) = [mean(x02_pf(1,:)) mean(x02_pf(2,:)) mean(x02_pf(3,:)) mean(x02_pf(4,:)) ...
                     mean(x02_pf(5,:)) mean(x02_pf(6,:)) mean(x02_pf(7,:)) mean(x02_pf(8,:)) ...
                     mean(x02_pf(9,:))]';
    x3_IMM(:,t+1) = [mean(x03_pf(1,:)) mean(x03_pf(2,:)) mean(x03_pf(3,:)) mean(x03_pf(4,:)) ...
                     mean(x03_pf(5,:)) mean(x03_pf(6,:)) mean(x03_pf(7,:)) mean(x03_pf(8,:)) ...
                     mean(x03_pf(9,:))]';
    
    S1 = R;
    S2 = R;
    S3 = R;
    for j = 1:N
       %求粒子滤波均值与各粒子值之差（9维）
       e1(:,j) = x1_IMM(:,t+1) - x01_pf(:,j);
       e2(:,j) = x2_IMM(:,t+1) - x02_pf(:,j);
       e3(:,j) = x3_IMM(:,t+1) - x03_pf(:,j);
       
       %更新各模型的状态量协方差矩阵（9*9） P = P0 + Pnew Pnew = e*e'
       P1_IMM(:,:,t+1) = P1_IMM(:,:,t+1)+e1(:,j)*e1(:,j)'*weight_model_1(j);
       P2_IMM(:,:,t+1) = P2_IMM(:,:,t+1)+e2(:,j)*e2(:,j)'*weight_model_2(j);
       P3_IMM(:,:,t+1) = P3_IMM(:,:,t+1)+e3(:,j)*e3(:,j)'*weight_model_3(j);
       
       %更新各模型观测值的协方差矩阵（3*3） S = R + Snew Snew = ero*ero'
       S1 = S1+ero1(:,j)*ero1(:,j)'*weight_model_1(j);
       S2 = S2+ero2(:,j)*ero2(:,j)'*weight_model_2(j);
       S3 = S3+ero3(:,j)*ero3(:,j)'*weight_model_3(j);
       
    end
    P1_IMM(:,:,t+1) = P1_IMM(:,:,t+1)+P1_IMM(:,:,t);
    P2_IMM(:,:,t+1) = P2_IMM(:,:,t+1)+P2_IMM(:,:,t);
    P3_IMM(:,:,t+1) = P3_IMM(:,:,t+1)+P3_IMM(:,:,t);
    
    %求各模型粒子滤波后平均的粒子滤波观测值与测量值之差
    r1_IMM = z(:,t+1)-H*x1_IMM(:,t+1);
    r2_IMM = z(:,t+1)-H*x2_IMM(:,t+1);
    r3_IMM = z(:,t+1)-H*x3_IMM(:,t+1);
                 

    %第三步 更新模型概率
    [u_IMM_N(:,t+1)] = Model_P_up(r1_IMM,r2_IMM,r3_IMM,S1,S2,S3,c_j);
    
    %第四步 模型输出综合
    [x_pro_IMMPF(:,t+1),P_IMMPF(:,:,t+1)] = Model_mix(x1_IMM(:,t+1),x2_IMM(:,t+1),x3_IMM(:,t+1),...
                                            P1_IMM(:,:,t+1),P2_IMM(:,:,t+1),P3_IMM(:,:,t+1),u_IMM_N(:,t+1));
    
    
    
end
toc

%==========================================================================
%计算均方根误差
%==========================================================================
z_IMMPF_sim = H*x_pro_IMMPF(:,:); 
Rmse_m = Compute_Rmse_z(z,z_true,simTime);%计算测量值的均方根误差
Rmse_IMMPF_sim = Compute_Rmse_z(z_IMMPF_sim,z_true,simTime);%计算测量值的均方根误差

%==========================================================================
%绘图
%==========================================================================
%模型概率
figure(1)
t  = 1:simTime;
plot(t,u_IMM_N(1,t),'k:',t,u_IMM_N(2,t),'r-',t,u_IMM_N(3,t),'b--','linewidth',2);grid on
title('IMM算法模型概率曲线');
xlabel('t/s');ylabel('模型概率');
legend('模型1','模型2','模型3');

figure(2)
plot3(z_true(1,:),z_true(2,:),z_true(3,:),'black');grid on;hold on;
plot3(z(1,:),z(2,:),z(3,:),'r');
plot3(z_IMMPF_sim(1,:),z_IMMPF_sim(2,:),z_IMMPF_sim(3,:),'b');
legend('真实值','测量值','IMMPF滤波值');
hold off
magnify
