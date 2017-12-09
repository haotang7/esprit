tic
clear variables;
close all;
j=sqrt(-1);
% ======= elevation and azimuth in degree ======== %
azimuth=[40;68];      %方位角[40;169];[42;149],[60;128] 
elevation=[30;55];   %俯仰角 [30;55];[55;78],[49;70]  
degrad=pi/180;%角度-弧度换算
azimuth1=azimuth*degrad;
elevation1=elevation*degrad;
% ===================signal and spherical array============================= %
snapshot=200;%快拍数
N_x=(1:snapshot);%采样点数； ·
I=size( azimuth, 1 );% number of   signal s
load mic_position;
theta_y=mic_position(:,1)*degrad;
fai_y=mic_position(:,2)*degrad;
M=32;
N=4;
N1=10;
kr=[4;4];
ka=kr;
L1=(N+1).^2;%Y矩阵的列长度
%%%%%%%%---------------------求出A矩阵--------------------------------%%%%%%%
y_omega=zeros(M,L1);%生成零矩阵Y(omega)
y_fai=zeros(1,L1);%生成零矩阵Y(fai)
for n1=0:N
    for m1=-n1:n1 
        for  k11=1:M
             y_omega(k11,n1.^2+n1+m1+1)=spherical_har(m1,n1,theta_y(k11),fai_y(k11));     % 生成Y(omega)矩阵
        end
        for i1=1:I                
            y_fai(i1,n1.^2+n1+m1+1)=spherical_har(m1,n1,elevation1(i1),azimuth1(i1));     % 生成Y(fai)矩阵
        end
    end
end 
bn=[];
        for n=0:N
            bn_temp=bn_rigid(kr,ka,n)*ones(1,2*n+1);
            bn=[bn bn_temp];    
        end
Bn=diag(bn(1,:));%生成Bn矩阵
Y_f=y_fai';
Y_fai=Bn*y_fai';%求Y（fai）的复共轭
A=y_omega* Y_fai;
[n,p]=size(A);
K=5;
l=25/K;
ss=randn(I,snapshot)+j*randn(I,snapshot);
snr=15;
refp=10.^(snr/20);
S=refp.*ss;
[p,T]=size(S);
Signal=A*S;
noise=randn(M,snapshot)+j*randn(M,snapshot);
x=Signal+noise;
C=pinv(y_omega);%tubu
C2=pinv(Bn);
X1=C2*C*x;
A1=C2*C*A;
W_ex=orth(A1);
W = zeros(l,p,K);
U = zeros(l,p,K);
Y_1 = zeros(2,T+N,K);
Z = zeros(p,p,K);
S = zeros(p,p,K);
P = zeros(l,l,K);
M = zeros(l,l,K);
x_ = zeros(l,N,K);
X =zeros(l,T,K);
X_=zeros(l,T+N,K);
y=zeros(p,1,K);
y1=zeros(p,T+N,K);
x=zeros(l,2,K);
Y=zeros(p,2,K);
H=zeros(p,2,K);
Gamma=zeros(2,2,K);
G=zeros(p,2,K);
E=zeros(l,2,K);
F=zeros(p,4,K);
Tmp=zeros(p,4,K);
Sigma=zeros(p,4,K);
Q=zeros(p,p,K);
O=zeros(4,4,K);
L=zeros(4,p,K);
T_=zeros(p,p,K);
P_=zeros(l,p,K);
J=zeros(4,4,K);
X11=zeros(l,T);
X_11=zeros(l,T+N);
W1 = zeros(l,p);
P1_=zeros(l,p,K);
T1=zeros(p,p,K);
for t=1:T
    if t==1
      parfor j=1:K
            X(:,:,j) = X1((j-1)*l+1:j * l,:);
            W(:,:,j)=[eye(p),zeros(p,l-p)].';
            U(:,:,j)=W(:,:,j);
            Z(:,:,j)=eye(p);
            S(:,:,j)=Z(:,:,j);
            P(:,:,j)=eye(l);
            M(:,:,j)=eye(l);
            Z_tmp=zeros(l-p,N-p);
            x_(:,:,j)=blkdiag(eye(p), Z_tmp);
            Y_1(:,:,j)=[eye(p),zeros(p,T+N-p)];
            X_(:,:,j)=[x_(:,:,j),X(:,:,j)];
       end
    end
    for i=1:K
        y(:,:,i)= W(:,:,i)' * X(:,t,i);
        y1(:,:,i)=[zeros(p,N+t-1),y(:,:,i),zeros(p,T-t)];
        Y_1(:,:,i) = Y_1(:,:,i)+y1(:,:,i);
        for k=1:t+N
            Y_11(:,k,i)=Y_1(:,k,i);
        end
        x(:,:,i) = [X(:,t,i) X_(:,t,i)];
        Y(:,:,i) = [y(:,:,i) Y_11(:,t,i)];
        H(:,:,i)=Z(:,:,i)*Y(:,:,i);
        Gamma(:,:,i)=(diag([1,-1])+H(:,:,i)'*Y(:,:,i))^(-1);
        G(:,:,i)=H(:,:,i)*Gamma(:,:,i);
        Z(:,:,i)=Z(:,:,i)-G(:,:,i)*H(:,:,i)';
        E(:,:,i)=x(:,:,i)-U(:,:,i)*Y(:,:,i);
        U(:,:,i)=U(:,:,i)+E(:,:,i)*G(:,:,i)';
        F(:,:,i)=[G(:,:,i) U(:,:,i)'*E(:,:,i)];
        Tmp(:,:,i)=S(:,:,i)'*F(:,:,i);
        [Q(:,:,i),Sigma(:,:,i),O(:,:,i)]=svd(Tmp(:,:,i));
        L(:,:,i)= O(:,1:2,i)*Sigma(:,1:2,i);
        p1 = size(Sigma(:,:,i), 1);
        J(:,:,i) = [E(:,:,i)'*E(:,:,i) eye(2); eye(2) zeros(2,2)];
        P_(:,:,i) = [E(:,:,i) zeros(l, 2)] * L(:,:,i);
        T_(:,:,i) = (eye(p1) + L(:,:,i)'*J(:,:,i)*L(:,:,i))^(-0.5) - eye(p1);
        S(:,:,i) = S(:,:,i) + S(:,:,i) * Q(:,:,i) * T_(:,:,i) * Q(:,:,i)';
        P1_(:,:,i) = W(:,:,i) * Q(:,:,i) * T_(:,:,i) + P_(:,:,i) * (eye(p1) + T_(:,:,i));
        W(:,:,i) = W(:,:,i) + P1_(:,:,i) * Q(:,:,i)';
        if i==1
           X11=X(:,:,i);
           X_11=X_(:,:,i);
           W1=W(:,:,i);
        end
        P(:,:,i)=P(:,:,i)+X(:,t,i)*X11(:,t)'-X_(:,t,i)*X_11(:,t)';
        M(:,:,i)=M(:,:,i)+X11(:,t)*X11(:,t)'-X_11(:,t)*X_11(:,t)';
        T1(:,:,i)=pinv(W(:,:,i))*(P(:,:,i)*M(:,:,i)^-1)*W1;
        if i==1
            W_1=W(:,:,i).';
        else
            W_1=[W_1,(W(:,:,i)*T1(:,:,i)).'];
        end
    end
    W_1=W_1.';
    I=eye(25);
    SEP(t)=trace(W_1'*(I-W_ex*W_ex')*W_1)/trace(W_1'*(W_ex*W_ex')*W_1);
end
toc
t=1:20:200;
semilogy(t,SEP(t),'m:s');
grid minor;