function [u] = Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
%ģ�͸��ʸ��º���

% u--ģ�͸���
% r1--ģ��1Ԥ�����
% r2--ģ��2Ԥ�����
% r3--ģ��3Ԥ�����
% S1--ģ��1Ԥ�����Э�������
% S2--ģ��2Ԥ�����Э�������
% S3--ģ��3Ԥ�����Э�������
% c_j--ģ�ͻ�ϸ���

%������Ȼ����
Lfun1 = (1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));
Lfun2 = (1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));
Lfun3 = (1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));

%��һ��
Lfun11 = Lfun1^2/(Lfun1+Lfun2+Lfun3);
Lfun21 = Lfun2^2/(Lfun1+Lfun2+Lfun3);
Lfun31 = Lfun3^2/(Lfun1+Lfun2+Lfun3);

%����ģ�͸��¸���
c = [Lfun11,Lfun21,Lfun31]*c_j;
%��һ��
u = (1/c).*[Lfun11,Lfun21,Lfun31]'.*c_j;