function [r,m]=mulGSND(F,x0,h,eps)

% F 目标函数
% x0 初始值
% h 
% eps 误差
format long;
if nargin==3
    eps=1.0e-8;
end
        syms fai x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25;
        Us=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12;x13;x14;x15;x16;x17;x18;x19;x20;x21;x22;x23;x24;x25];
        v=[fai;Us];
n = length(x0);
% x0 = transpose(x0);
m=1;
tol=1;
while tol>eps
    fx = subs(F,findsym(F),x0);
    J = zeros(n,n);
%     J=jacobian(F,v); 
    for i=1:n
        x1 = x0;
        x1(i) = x1(i)+h;
        J(:,i) = (subs(F,findsym(F),x1)-fx)/h;
    end
    DF = pinv(transpose(J)*J)*transpose(J);
    r=x0-DF*fx;  %核心迭代公式
    tol=norm(r-x0);
    x0 = r;
    m=m+1;
    if(m>100)                                              %迭代步数控制
        disp('迭代步数太多，可能不收敛！');
        return;
    end
end
format short;