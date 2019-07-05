function [w,x_pf,e] = Particle_Residual_sim(N,x,z,F,H,r,q,d)
%�����˲� IMMPF_simר�ð�
R = r*eye(3);
x_pre = zeros(size(x));
x_pf = zeros(size(x));
weight = zeros(1,N);
z_pre = zeros(size(z));
randn('state',sum(100*clock));
R1 = sqrt(r)*randn(3,N);
Q1 = sqrt(q)*randn(9,N);
e = zeros(3,N); %�����˲�Ԥ��ֵ�����ֵ��Ĳв�
%����ÿ�����ӵ�Ȩ��
for i = 1:N
    x_pre(:,i) = F*x(:,i)+Q1(:,i);
    z_pre(:,i) = H*x_pre(:,i);
    weight(:,i) = exp(-0.5*((z_pre(:,i)-z(:,:)))'*inv(R)*((z_pre(:,i)-z(:,:))))+0.1;
end
%��һ��
sum_weight = sum(weight);
w = zeros(size(weight));
for k = 1:N
    w(:,k) = weight(:,k)./sum_weight;
    
end
%�ز��� ����ʽ�ز���
outIndex = multinomialR(w(1,:)');
for j = 1:N
    x_pf(:,j) = x_pre(:,outIndex(1,j));
end

for i = 1:N
   e(:,i) = z -  H*x_pf(:,i);
end

end