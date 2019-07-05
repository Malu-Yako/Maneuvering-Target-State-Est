%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��Ȩ������
%     ���������ϸ����ע����ο�
%     ��Сƽ�����ң�������.�����˲�ԭ��Ӧ��[M].���ӹ�ҵ�����磬2017.4
%     ������ԭ�����+����+����+����ע��
%     ����˳����д��������ʾ�޸�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ʽ�ز����Ӻ���
% ���������weightΪԭʼ���ݶ�Ӧ��Ȩ�ش�С
% ���������outIndex�Ǹ���weightɸѡ�͸��ƽ��
function outIndex = multinomialR(weight);
Col=length(weight)
N_babies= zeros(1,Col);

 
cdf= cumsum(weight);
 
u=rand(1,Col)

 
uu=u.^(1./(Col:-1:1))
 
ArrayTemp=cumprod(uu)
 
u = fliplr(ArrayTemp);
j=1;
for i=1:Col
 
    while (u(i)>cdf(j))
        j=j+1;
    end
    N_babies(j)=N_babies(j)+1;
end;
index=1;
for i=1:Col
    if (N_babies(i)>0)
        for j=index:index+N_babies(i)-1
            outIndex(j) = i;
        end;
    end;
    index= index+N_babies(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
