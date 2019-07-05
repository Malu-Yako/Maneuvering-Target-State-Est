function [X,P,e,S] = Kalman(X_Forward,P_Forward,Z,F,G,Q,H,R)
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

%预测

X_Pre = F*X_Forward;
P_Pre = F*P_Forward*F'+G*Q*G';
%P_Pre = F*P_Forward*F'+Q;

%卡尔曼增益
e = Z-H*(F*X_Forward);
S = H*P_Pre*H'+R;

K = P_Pre*H'*inv(S)';

%修正滤波值与误差协方差矩阵
X = X_Pre+K*e;
M = K*H;
n = size(M);
I = eye(n);
P = (I-K*H)*P_Pre*(I-K*H)'+K*R*K';







end

