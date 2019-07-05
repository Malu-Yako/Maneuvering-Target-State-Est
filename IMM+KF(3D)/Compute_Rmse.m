function [Rmse] = Compute_Rmse(X,Z,simTime)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
R_1 = zeros(1,simTime);
R_2 = zeros(1,simTime);
R_3 = zeros(1,simTime);
for t = 1:simTime
    R_1(1,t) = (Z(1,t)-X(1,t))^2;
    R_2(1,t) = (Z(2,t)-X(4,t))^2;
    R_3(1,t) = (Z(3,t)-X(7,t))^2;
end
Rmse_1 = (mean(R_1))^0.5;
Rmse_2 = (mean(R_2))^0.5;
Rmse_3 = (mean(R_3))^0.5;

Rmse = (Rmse_1+Rmse_2+Rmse_3)/3;
end

