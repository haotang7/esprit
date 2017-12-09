function [r,n]=mulNewton(F,x0,eps)
%非线性方程组：F
%初始迭代值：x0
%解的精度：eps
%求得的一组解：r
%迭代步数：n

if nargin==2
    eps=1.0e-6;
end
x0=transpose(x0);
Fx=subs(F,findsym(F),x0);    %F函数就是S目标函数
dF=jacobian(F);              %对F求Jacobian矩阵
dFx=subs(dF,findsym(dF),x0); 
r=x0-Fx\dFx;                 %反映高斯迭代的核心迭代公式
n=1;
to1=1;
while to1>eps
% while n>5
    x0=r;
    Fx=subs(F,findsym(F),x0);
    dFx=subs(dF,findsym(dF),x0);
    r=x0-Fx\dFx;              %核心迭代公式
    to1=norm(r-x0);
    n=n+1;
    if(n>100000)
        disp('迭代次数太多，可能不收敛！');
        return;
    end
end    