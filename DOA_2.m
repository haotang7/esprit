clear variables;
j=sqrt(-1);
% ======= elevation and azimuth in degree ======== %
azimuth=[42;149];      %方位角[40;169];[42;149],[60;128] 
elevation=[55;78];   %俯仰角 [30;55];[55;78],[49;70]  
degrad=pi/180;%角度-弧度换算
azimuth1=azimuth*degrad;
elevation1=elevation*degrad;
% ===================signal and spherical array============================= %
snapshot=2000;%快拍数
N_x=(1:snapshot);%采样点数； ·
I=size( azimuth, 1 );% number of   signal s
load mic_position;
theta_y=mic_position(:,1)*degrad;
fai_y=mic_position(:,2)*degrad;
M=32;
N=4;
N1=200;
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
K=5;
l=25/K;
ss=randn(I,snapshot);
monte_num=100;
snr=15;
refp=10.^(snr/20);
S=refp.*ss;
[q,T]=size(S);
Signal=A*S;
azim_estimated_c=zeros(q,snapshot,100);
elev_estimated_c=zeros(q,snapshot,100);
for p=1:100
    try
        noise=randn(M,snapshot);
        x_y=Signal+noise;
        C=pinv(y_omega);
        C2=pinv(Bn);
        X1=C2*C*x_y;
        W = zeros(l,q,K);
        U = zeros(l,q,K);
        Y_1 = zeros(2,T+N1,K);
        Z = zeros(q,q,K);
        S = zeros(q,q,K);
        P = zeros(l,l,K);
        M1 = zeros(l,l,K);
        x_ = zeros(l,N1,K);
        X =zeros(l,T,K);
        X_=zeros(l,T+N1,K);
        y=zeros(q,1,K);
        y1=zeros(q,T+N1,K);
        x=zeros(l,2,K);
        Y=zeros(q,2,K);
        H=zeros(q,2,K);
        Gamma=zeros(2,2,K);
        G=zeros(q,2,K);
        E=zeros(l,2,K);
        F=zeros(q,4,K);
        Tmp=zeros(q,4,K);
        Sigma=zeros(q,4,K);
        Q=zeros(q,q,K);
        O=zeros(4,4,K);
        L=zeros(4,q,K);
        T_=zeros(q,q,K);
        P_=zeros(l,q,K);
        J=zeros(4,4,K);
        X11=zeros(l,T);
        X_11=zeros(l,T+N1);
        W1 = zeros(l,q);
        P1_=zeros(l,q,K);
        for j=1:K
            X(:,:,j) = X1((j-1)*l+1:j * l,:);
            W(:,:,j)=[eye(q),zeros(q,l-q)].';
            U(:,:,j)=W(:,:,j);
            Z(:,:,j)=eye(q);
            S(:,:,j)=Z(:,:,j);
            P(:,:,j)=eye(l);
            M1(:,:,j)=eye(l);
            Z_tmp=zeros(l-q,N1-q);
            x_(:,:,j)=blkdiag(eye(q), Z_tmp);
            Y_1(:,:,j)=[eye(q),zeros(q,T+N1-q)];
            X_(:,:,j)=[x_(:,:,j),X(:,:,j)];
        end
        for t=1:snapshot
            for i=1:K
                y(:,:,i)= W(:,:,i)' * X(:,t,i);
                y1(:,:,i)=[zeros(q,N1+t-1),y(:,:,i),zeros(q,T-t)];
                Y_1(:,:,i) = Y_1(:,:,i)+y1(:,:,i);
                Y_11=zeros(q,t+N1,K);
                Y_11(:,:,i)=Y_1(:,1:t+N1,i);
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
                    X_11=[x_(:,:,i),X(:,:,i)];
                    W1=W(:,:,i);
                end
                P(:,:,i)=P(:,:,i)+X(:,t,i)*X11(:,t)'-X_(:,t,i)*X_11(:,t)';
                M1(:,:,i)=M1(:,:,i)+X11(:,t)*X11(:,t)'-X_11(:,t)*X_11(:,t)';
                T1=pinv(W(:,:,i))*(P(:,:,i)*M1(:,:,i)^-1)*W1;
                if i==1
                    W_1=W(:,:,i).';
                else
                    W_1=[W_1,(W(:,:,i)*T1).'];
                end
            end
            W_1=W_1.';
            Sspace=W_1;
        %%%%%%%%%%%%%%%%%%%%%%%%% complex spherical harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            S1_c=Sspace([3,6,7,8,11,12,13,14,15],:);
            S0_c=Sspace([1,2,3,4,5,6,7,8,9],:);
            S2_c=Sspace([7,12,13,14,19,20,21,22,23],:);
            lamda1_c=[lamda(1,0),lamda(2,-1),lamda(2,0),lamda(2,1),lamda(3,-2),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(3,2)];
            lamda2_c=[lamda(2,0),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(4,-2),lamda(4,-1),lamda(4,0),lamda(4,1),lamda(4,2)];
            La_c=diag(lamda1_c);
            La2_c=diag(lamda2_c);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            S_a0_c=Sspace([4,8,9,12,13,14,15,16],:);
            S_a1_c=Sspace([1,3,4,5,6,7,8,9],:);
            S_a2_c=Sspace([7,13,14,19,20,21,22,23],:);
       
            mu1_c=[mu(1,1),mu(2,1),mu(2,2),mu(3,-1),mu(3,0),mu(3,1),mu(3,2),mu(3,3)];
            mu2_c=[mu(2,0),mu(3,0),mu(3,-1),mu(4,2),mu(4,1),mu(4,0),mu(4,-1),mu(4,-2)];          
            Mu1_c=diag(mu1_c);
            Mu2_c=diag(mu2_c);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT  %%%%%%%%%%%%%%%%%%%%%%%% 
            SS1_c=pinv(S1_c'*S1_c);
            FAI_c=SS1_c*S1_c'*(La_c*S0_c+La2_c*S2_c);
            [Vm_c,Dm_c]=eig(FAI_c);

            elev_estimated_c(:,t,p)=acos(diag(Dm_c))*180/pi;

      %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%%  
            SS_a11_c=pinv(S_a0_c'*S_a0_c);
            FAI_a1_c=SS_a11_c*S_a0_c'*(-Mu1_c*S_a1_c+Mu2_c*S_a2_c);
            [Vm_c1,Dm_c1]=eig(FAI_a1_c);

            azim_estimated_c(:,t,p)=angle(diag(Dm_c1)./(sin(acos(diag(Dm_c)))))*180/pi;  

        end
    catch
        p=p-1;
    end
end 

for t=1:snapshot
    for p=1:100
        elev_estimated_c_1(:,p)=elev_estimated_c(:,t,p);
        elev_estimated_c_1(:,p)=sort(elev_estimated_c_1(:,p));
        azim_estimated_c_1(:,p)=azim_estimated_c(:,t,p);
        azim_estimated_c_1(:,p)=sort(azim_estimated_c_1(:,p));
    end
    RMSE_azim_esprit_c(:,t)=sqrt(sum((((azim_estimated_c_1-azimuth*ones(1,p))).^2)')/p);
    RMSE_elev_estimated_c(:,t)=sqrt(sum((abs(((elev_estimated_c_1-elevation*ones(1,p)))).^2)')/p);
end

RMSE_azim_E_c=sum(RMSE_azim_esprit_c)/2;
RMSE_elev_E_c=sum(RMSE_elev_estimated_c)/2;

t=1:T;

figure;
plot(t,RMSE_elev_E_c,'b-^');
legend('D-SHESPRIT');
xlabel('t');
ylabel('RMSE of elevation (degree)');

figure;
plot(t,RMSE_azim_E_c,'b-^');
legend('D-SHESPRIT');
xlabel('t');
ylabel('RMSE of azimuth (degree)');
