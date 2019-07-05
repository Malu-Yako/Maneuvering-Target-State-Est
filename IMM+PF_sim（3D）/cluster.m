function [index_cluster,Particle_X] = cluster(Particle_X,x,r,N)
%���� ��������� ����һʱ�������˲�������N�����ӣ��������о���Ĳ���ֵ������ľ��룬������
%   �˴���ʾ��ϸ˵��
rr = zeros(9,N);
R = zeros(1,N);
index_cluster = zeros(1,N);
V = [0.2 0.06 0.06 0.2 0.06 0.06 0.2 0.06 0.06];%����Ȩֵ
%����ÿ�����������ֵ�Ĳ�ֵ�ľ���ֵ
for i = 1:N
    rr(:,i) = abs(Particle_X(:,i)-x);
    R(i) = V*rr(:,i);
end
%�������ϱ�� ��ֵ�����������ı�־1 �ò���ֵ�����滻
for i = 1:N
    if R(i)>=r
        index_cluster(i) = 1;
        Particle_X(:,i) = x;
    end
end



