%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 版权声明：
%     本程序的详细中文注释请参考
%     黄小平，王岩，缪鹏程.粒子滤波原理及应用[M].电子工业出版社，2017.4
%     书中有原理介绍+例子+程序+中文注释
%     如果此程序有错误，请对提示修改
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 多项式重采样子函数
% 输入参数：weight为原始数据对应的权重大小
% 输出参数：outIndex是根据weight筛选和复制结果
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
