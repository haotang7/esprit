tic
clear variables;
close all;
p=2;
n=48000;
K=4;
T=3000;
q=4;
N=500;
l=n/K;
SNR=15;
A=randn(n,p);
S=randn(p,T);
mapminmax(S,0,1);
X_1=A*S;  
X = awgn(X_1,SNR,'measured');
A1=orth(A);  
       W=[eye(p),zeros(p,n-p)].';
       U=W;
       Z=eye(p);
       S=Z;
       Z_tmp = zeros(n-p,N-p);
       x_=blkdiag(eye(p), Z_tmp);
       Y_=[eye(p),zeros(p,N-p)];
       X_ = [x_ X];
for t=1:T
       y = W' * X(:, t);
       Y_ = [Y_ y];
       x = [X(:,t) X_(:,t)];
       Y = [y Y_(:,t)];
      H=Z*Y;
    Gamma=(diag([1,-1])+H'*Y)^(-1);
    G=H*Gamma;
    Z=Z-G*H';
    E=x-U*Y;
    U=U+E*G';
    F=[G U'*E];
    Tmp=S'*F;
      [Q,Sigma,O]=svd(Tmp);
    L= O(:,1:2)*Sigma(:,1:2);
    p1 = size(Sigma, 1);
    J = [E' * E eye(2); eye(2) zeros(2,2)];
    P_ = [E zeros(n, 2)] * L;
    T = (eye(p1) + L'*J*L)^(-0.5) - eye(p1);
    S = S + S * Q * T * Q';
    P1_ = W * Q * T + P_ * (eye(p1) + T);
    W = W + P1_ * Q';
end
toc
