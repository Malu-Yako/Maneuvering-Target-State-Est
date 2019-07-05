function [index_cluster,Particle_X] = cluster(Particle_X,x,r,N)
%聚类 按距离聚类 从上一时刻粒子滤波采样的N个粒子，用来进行聚类的测量值，容许的距离，粒子数
%   此处显示详细说明
rr = zeros(9,N);
R = zeros(1,N);
index_cluster = zeros(1,N);
V = [0.2 0.06 0.06 0.2 0.06 0.06 0.2 0.06 0.06];%误差的权值
%计算每个粒子与测量值的差值的绝对值
for i = 1:N
    rr(:,i) = abs(Particle_X(:,i)-x);
    R(i) = V*rr(:,i);
end
%给粒子上标记 差值超过容许距离的标志1 用测量值进行替换
for i = 1:N
    if R(i)>=r
        index_cluster(i) = 1;
        Particle_X(:,i) = x;
    end
end



