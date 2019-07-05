function [X,P,e,S] = U_Kalman(X_Forward,P_Forward,Z,F,G,Q,H,R)
%UNTITLED2 此处显示有关此函数的摘要
%Z--观测数据矢量
%F--系统模型状态转移矩阵
%G--系统模型噪声系数矩阵
%Q--系统模型噪声方差
%H--系统模型量测矩阵
%R--量测模型协方差矩阵
%X_Forward--上一时刻的值
%P_Forward

%X--输出估计状态矢量
%P--输出估计状态协方差矩阵
%e--残差
%S--残差协方差矩阵

%==========================================================================
%U_Kalman
%==========================================================================

%UT变换
L = 9; %状态维数
alpha = 1;
kalpha = 0;
belta = 2;
ramda = 3-L;
for j = 1:2*L+1
    Wm(j) = 1/(2*(L+ramda));
    Wc(j) = 1/(2*(L+ramda));
end
Wm(1) = ramda/(L+ramda);%求期望时的权重
Wc(1) = ramda/(L+ramda)+1-alpha^2+belta;

%第一步 获得一组sigma点集
cho = (chol(P_Forward*(L+ramda)))';
for k = 1:L
    xgama_P1(:,k) = X_Forward+cho(:,k);
    xgama_P2(:,k) = X_Forward-cho(:,k);
end
X_Sigma = [X_Forward,xgama_P1,xgama_P2] %Sigma点集

%第二步 对Sigma点集进行一步预测
X_Sigma_pre = F*X_Sigma;

%第三步 利用第二步的结果计算均值与协方差
X_pre = zeros(4,1); %预测均值
for k = 1:2*L+1
    X_pre = X_pre+Wm(k)*X_Sigma_pre(:,k);
end
P_pre = zeros(4,4);%预测协方差
for k = 1:2*L+1
    P_pre = P_pre + Wc(k)*(X_Sigma_pre(:,k)-X_pre)*(X_Sigma_pre(:,k)-X_pre)';
end
P_pre = P_pre +G*Q*G';
%P_pre = P_pre +Q;
%第四步：根据预测值，再进行依次UT变换 得到新的sigma点集
chol2 = (chol((L+ramda)*P_pre))';
for k = 1:L
    xgama2_P1(:,k) = X_pre+chol2(:,k);
    xgama2_P2(:,k) = X_pre-chol2(:,k);
end
X2_Sigma = [X_pre xgama2_P1 xgama2_P2];

%第五步：观测预测
Z_Sigma_pre = zeros(2,9);
for k = 1:2*L+1 %预测的观测值
    Z_Sigma_pre(:,k)=H*X2_Sigma(:,k);
end

%第六步：计算预测的观测值的均值与协方差
Z_pre = 0;
for k = 1:2*L+1
    Z_pre = Z_pre+Wm(k)*Z_Sigma_pre(:,k);
end
Pzz = zeros(2,2); %观测值的协方差矩阵
for k = 1:2*L+1
    Pzz = Pzz+Wc(k)*(Z_Sigma_pre(:,k)-Z_pre)*(Z_Sigma_pre(:,k)-Z_pre)';
end
Pzz = Pzz+R;
    
Pxz = 0;%状态值的协方差矩阵
for k = 1:2*L+1
    Pxz = Pxz+Wc(k)*(X2_Sigma(:,k)-X_pre)*(Z_Sigma_pre(:,k)-Z_pre)';
end

%第七步：计算卡尔曼滤波增益
K = Pxz*inv(Pzz);

%第八步：状态 方差 更新
X = X_pre+K*(Z-Z_pre);
P = P_pre-K*Pzz*K';
S = H*P_pre*H'+R;
e = Z-H*(F*X_Forward);
    






















