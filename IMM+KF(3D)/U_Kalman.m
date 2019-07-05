function [X,P,e,S] = U_Kalman(X_Forward,P_Forward,Z,F,G,Q,H,R)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%Z--�۲�����ʸ��
%F--ϵͳģ��״̬ת�ƾ���
%G--ϵͳģ������ϵ������
%Q--ϵͳģ����������
%H--ϵͳģ���������
%R--����ģ��Э�������
%X_Forward--��һʱ�̵�ֵ
%P_Forward

%X--�������״̬ʸ��
%P--�������״̬Э�������
%e--�в�
%S--�в�Э�������

%==========================================================================
%U_Kalman
%==========================================================================

%UT�任
L = 9; %״̬ά��
alpha = 1;
kalpha = 0;
belta = 2;
ramda = 3-L;
for j = 1:2*L+1
    Wm(j) = 1/(2*(L+ramda));
    Wc(j) = 1/(2*(L+ramda));
end
Wm(1) = ramda/(L+ramda);%������ʱ��Ȩ��
Wc(1) = ramda/(L+ramda)+1-alpha^2+belta;

%��һ�� ���һ��sigma�㼯
cho = (chol(P_Forward*(L+ramda)))';
for k = 1:L
    xgama_P1(:,k) = X_Forward+cho(:,k);
    xgama_P2(:,k) = X_Forward-cho(:,k);
end
X_Sigma = [X_Forward,xgama_P1,xgama_P2] %Sigma�㼯

%�ڶ��� ��Sigma�㼯����һ��Ԥ��
X_Sigma_pre = F*X_Sigma;

%������ ���õڶ����Ľ�������ֵ��Э����
X_pre = zeros(4,1); %Ԥ���ֵ
for k = 1:2*L+1
    X_pre = X_pre+Wm(k)*X_Sigma_pre(:,k);
end
P_pre = zeros(4,4);%Ԥ��Э����
for k = 1:2*L+1
    P_pre = P_pre + Wc(k)*(X_Sigma_pre(:,k)-X_pre)*(X_Sigma_pre(:,k)-X_pre)';
end
P_pre = P_pre +G*Q*G';
%P_pre = P_pre +Q;
%���Ĳ�������Ԥ��ֵ���ٽ�������UT�任 �õ��µ�sigma�㼯
chol2 = (chol((L+ramda)*P_pre))';
for k = 1:L
    xgama2_P1(:,k) = X_pre+chol2(:,k);
    xgama2_P2(:,k) = X_pre-chol2(:,k);
end
X2_Sigma = [X_pre xgama2_P1 xgama2_P2];

%���岽���۲�Ԥ��
Z_Sigma_pre = zeros(2,9);
for k = 1:2*L+1 %Ԥ��Ĺ۲�ֵ
    Z_Sigma_pre(:,k)=H*X2_Sigma(:,k);
end

%������������Ԥ��Ĺ۲�ֵ�ľ�ֵ��Э����
Z_pre = 0;
for k = 1:2*L+1
    Z_pre = Z_pre+Wm(k)*Z_Sigma_pre(:,k);
end
Pzz = zeros(2,2); %�۲�ֵ��Э�������
for k = 1:2*L+1
    Pzz = Pzz+Wc(k)*(Z_Sigma_pre(:,k)-Z_pre)*(Z_Sigma_pre(:,k)-Z_pre)';
end
Pzz = Pzz+R;
    
Pxz = 0;%״ֵ̬��Э�������
for k = 1:2*L+1
    Pxz = Pxz+Wc(k)*(X2_Sigma(:,k)-X_pre)*(Z_Sigma_pre(:,k)-Z_pre)';
end

%���߲������㿨�����˲�����
K = Pxz*inv(Pzz);

%�ڰ˲���״̬ ���� ����
X = X_pre+K*(Z-Z_pre);
P = P_pre-K*Pzz*K';
S = H*P_pre*H'+R;
e = Z-H*(F*X_Forward);
    






















