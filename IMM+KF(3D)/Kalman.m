function [X,P,e,S] = Kalman(X_Forward,P_Forward,Z,F,G,Q,H,R)
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

%Ԥ��

X_Pre = F*X_Forward;
P_Pre = F*P_Forward*F'+G*Q*G';
%P_Pre = F*P_Forward*F'+Q;

%����������
e = Z-H*(F*X_Forward);
S = H*P_Pre*H'+R;

K = P_Pre*H'*inv(S)';

%�����˲�ֵ�����Э�������
X = X_Pre+K*e;
M = K*H;
n = size(M);
I = eye(n);
P = (I-K*H)*P_Pre*(I-K*H)'+K*R*K';







end

