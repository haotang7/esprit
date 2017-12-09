function [b] = Trans3(x)

%   变换函数   导向向量仅仅包含球谐函数
%   x   输入维度
%   
%   b   输出的变换矩阵
b=zeros(x,x);
if(mod(x,2)==1)         %判断X是奇数还是偶数
b1=sqrt(2);
b2=Trans4(3);
b3=Trans4(5);
b4=Trans4(7);
b5=Trans4(9);
end
b=blkdiag(b1,b2,b3,b4,b5);
b=1./sqrt(2).*b;

function [a] = Trans4(k)
kk=fix(k/2);
j=sqrt(-1);
a=zeros(k,k);

a11=eye(kk);
a12=zeros(kk,1);
a13=flipud(Trans5(kk));
a21=zeros(1,kk);
a22=sqrt(2);
a23=zeros(1,kk);
a31=-j*inv_eye(kk);
a32=zeros(kk,1);
a33=j*Trans5(kk);
a=[a11 a12 a13;a21 a22 a23;a31 a32 a33];

function [cc] = Trans5(kk)
for i=1:kk
    c(i)=(-1)^i;
end
cc=diag(c);



    
    