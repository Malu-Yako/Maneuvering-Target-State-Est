function [u] = Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
%模型概率更新函数

% u--模型概率
% r1--模型1预测误差
% r2--模型2预测误差
% r3--模型3预测误差
% S1--模型1预测误差协方差矩阵
% S2--模型2预测误差协方差矩阵
% S3--模型3预测误差协方差矩阵
% c_j--模型混合概率

%计算似然函数
Lfun1 = (1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));
Lfun2 = (1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));
Lfun3 = (1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));

%归一化
Lfun11 = Lfun1^2/(Lfun1+Lfun2+Lfun3);
Lfun21 = Lfun2^2/(Lfun1+Lfun2+Lfun3);
Lfun31 = Lfun3^2/(Lfun1+Lfun2+Lfun3);

%计算模型更新概率
c = [Lfun11,Lfun21,Lfun31]*c_j;
%归一化
u = (1/c).*[Lfun11,Lfun21,Lfun31]'.*c_j;