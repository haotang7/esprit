clear variables;
tic
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
W_ex=orth(A);
[n,q]=size(A);
K=5;
l=25/K;
ss=randn(I,snapshot)+j*randn(I,snapshot);
monte_num=100;
snr=0:5:20;
startmatlabpool(4)
for ii=1:length(snr)
    refp=10.^(snr(ii)/20);
    S=refp.*ss;
    [q,T]=size(S);
    Signal=A*S;
    for p=1:monte_num
        noise=randn(M,snapshot)+j*randn(M,snapshot);
        X1=Signal+noise;
        C=pinv(y_omega);%tub
        C2=pinv(Bn);
        x_y=C2*C*X1;
        tic
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
        T1=zeros(q,q,K);
        for t=1:T
            if t==1
                parfor j=1:K
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
            end
            y(:,:,1)= W(:,:,1)' * X(:,t,1);
            y1(:,:,1)=[zeros(q,N1+t-1),y(:,:,1),zeros(q,T-t)];
            Y_1(:,:,1) = Y_1(:,:,1)+y1(:,:,1);
            for k=1:t+N1
                Y_11(:,k,1)=Y_1(:,k,1);
            end
            x(:,:,1) = [X(:,t,1) X_(:,t,1)];
            Y(:,:,1) = [y(:,:,1) Y_11(:,t,1)];
            H(:,:,1)=Z(:,:,1)*Y(:,:,1);
            Gamma(:,:,1)=(diag([1,-1])+H(:,:,1)'*Y(:,:,1))^(-1);
            G(:,:,1)=H(:,:,1)*Gamma(:,:,1);
            Z(:,:,1)=Z(:,:,1)-G(:,:,1)*H(:,:,1)';
            E(:,:,1)=x(:,:,1)-U(:,:,1)*Y(:,:,1);
            U(:,:,1)=U(:,:,1)+E(:,:,1)*G(:,:,1)';
            F(:,:,1)=[G(:,:,1) U(:,:,1)'*E(:,:,1)];
            Tmp(:,:,1)=S(:,:,1)'*F(:,:,1);
            [Q(:,:,1),Sigma(:,:,1),O(:,:,1)]=svd(Tmp(:,:,1));
            L(:,:,1)= O(:,1:2,1)*Sigma(:,1:2,1);
            p1 = size(Sigma(:,:,1), 1);
            J(:,:,1) = [E(:,:,1)'*E(:,:,1) eye(2); eye(2) zeros(2,2)];
            P_(:,:,1) = [E(:,:,1) zeros(l, 2)] * L(:,:,1);
            T_(:,:,1) = (eye(p1) + L(:,:,1)'*J(:,:,1)*L(:,:,1))^(-0.5) - eye(p1);
            S(:,:,1) = S(:,:,1) + S(:,:,1) * Q(:,:,1) * T_(:,:,1) * Q(:,:,1)';
            P1_(:,:,1) = W(:,:,1) * Q(:,:,1) * T_(:,:,1) + P_(:,:,1) * (eye(p1) + T_(:,:,1));
            W(:,:,1) = W(:,:,1) + P1_(:,:,1) * Q(:,:,1)';
            X11=X(:,:,1);
            X_11=X_(:,:,1);
            W1=W(:,:,1);
            W_1=W(:,:,1).';
            parfor i=2 :K
                y(:,:,i)= W(:,:,i)' * X(:,t,i);
                y1(:,:,i)=[zeros(q,N1+t-1),y(:,:,i),zeros(q,T-t)];
                Y_1(:,:,i) = Y_1(:,:,i)+y1(:,:,i);
            end
            parfor i=2 :K
                for k=1:t+N1
                    Y_11(:,k,i)=Y_1(:,k,i);
                end
            end
            parfor i=2 :K
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
            end
            parfor i=2 :K
                for k=1:2
                    O1(:,k,i)= O(:,k,i);
                    Sigma1(:,k,i)=Sigma(:,k,i);
                end
            end
            parfor i=2 :K
                L(:,:,i)= O1(:,:,i)*Sigma1(:,:,i);
                p1 = size(Sigma1(:,:,i), 1);
                J(:,:,i) = [E(:,:,i)'*E(:,:,i) eye(2); eye(2) zeros(2,2)];
                P_(:,:,i) = [E(:,:,i) zeros(l, 2)] * L(:,:,i);
                T_(:,:,i) = (eye(p1) + L(:,:,i)'*J(:,:,i)*L(:,:,i))^(-0.5) - eye(p1);
                S(:,:,i) = S(:,:,i) + S(:,:,i) * Q(:,:,i) * T_(:,:,i) * Q(:,:,i)';
                P1_(:,:,i) = W(:,:,i) * Q(:,:,i) * T_(:,:,i) + P_(:,:,i) * (eye(p1) + T_(:,:,i));
                W(:,:,i) = W(:,:,i) + P1_(:,:,i) * Q(:,:,i)';
                P(:,:,i)=P(:,:,i)+X(:,t,i)*X11(:,t)'-X_(:,t,i)*X_11(:,t)';
                M1(:,:,i)=M1(:,:,i)+X11(:,t)*X11(:,t)'-X_11(:,t)*X_11(:,t)';
                T1(:,:,i)=pinv(W(:,:,i))*(P(:,:,i)*M1(:,:,i)^-1)*W1;
                W_1=[W_1,(W(:,:,i)*T1(:,:,i)).'];
            end
            W_1=W_1.';
        end
        Sspace_GMNS=W_1;
%%%%%%%%%%%%%%%%%%%%%%%%% complex spherical harmonics GMNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        S1_c_GMNS=Sspace_GMNS([3,6,7,8,11,12,13,14,15],:);
        S0_c_GMNS=Sspace_GMNS([1,2,3,4,5,6,7,8,9],:);
        S2_c_GMNS=Sspace_GMNS([7,12,13,14,19,20,21,22,23],:);
        lamda1_c_GMNS=[lamda(1,0),lamda(2,-1),lamda(2,0),lamda(2,1),lamda(3,-2),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(3,2)];
        lamda2_c_GMNS=[lamda(2,0),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(4,-2),lamda(4,-1),lamda(4,0),lamda(4,1),lamda(4,2)];
        La_c_GMNS=diag(lamda1_c_GMNS);
        La2_c_GMNS=diag(lamda2_c_GMNS);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        S_a0_c_GMNS=Sspace_GMNS([4,8,9,12,13,14,15,16],:);
        S_a1_c_GMNS=Sspace_GMNS([1,3,4,5,6,7,8,9],:);
        S_a2_c_GMNS=Sspace_GMNS([7,13,14,19,20,21,22,23],:);
       
        mu1_c_GMNS=[mu(1,1),mu(2,1),mu(2,2),mu(3,-1),mu(3,0),mu(3,1),mu(3,2),mu(3,3)];
        mu2_c_GMNS=[mu(2,0),mu(3,0),mu(3,-1),mu(4,2),mu(4,1),mu(4,0),mu(4,-1),mu(4,-2)];          
        Mu1_c_GMNS=diag(mu1_c_GMNS);
        Mu2_c_GMNS=diag(mu2_c_GMNS);
        %%%%%%%%%%%%%%%%%%%%%%%%% complex spherical harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT  %%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT GMNS %%%%%%%%%%%%%%%%%%%%%%%% 
        SS1_c_GMNS=pinv(S1_c_GMNS'*S1_c_GMNS);
        FAI_c_GMNS=SS1_c_GMNS*S1_c_GMNS'*(La_c_GMNS*S0_c_GMNS+La2_c_GMNS*S2_c_GMNS);
        [Vm_c_GMNS,Dm_c_GMNS]=eig(FAI_c_GMNS);

        elev_estimated_c_GMNS(:,p)=acos(diag(Dm_c_GMNS))*180/pi;
        elev_estimated_c_GMNS=sort(elev_estimated_c_GMNS);
      %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%%   
       %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT GMNS %%%%%%%%%%%%%%%%%%%%%%%%%%      
        SS_a11_c_GMNS=pinv(S_a0_c_GMNS'*S_a0_c_GMNS);
        FAI_a1_c_GMNS=SS_a11_c_GMNS*S_a0_c_GMNS'*(-Mu1_c_GMNS*S_a1_c_GMNS+Mu2_c_GMNS*S_a2_c_GMNS);
        [Vm_c1,Dm_c1_GMNS]=eig(FAI_a1_c_GMNS);

        azim_estimated_c_GMNS(:,p)=angle(diag(Dm_c1_GMNS)./(sin(acos(diag(Dm_c_GMNS)))))*180/pi;  
        azim_estimated_c_GMNS=sort(azim_estimated_c_GMNS);      
    end
    RMSE_azim_esprit_c_GMNS(:,ii)=sqrt(sum((((azim_estimated_c_GMNS-azimuth*ones(1,p))).^2)')/p);
    RMSE_elev_estimated_c_GMNS(:,ii)=sqrt(sum((abs(((elev_estimated_c_GMNS-elevation*ones(1,p)))).^2)')/p);
     
end

RMSE_azim_E_c_GMNS=sum(RMSE_azim_esprit_c_GMNS)/2;
RMSE_elev_E_c_GMNS=sum(RMSE_elev_estimated_c_GMNS)/2;


figure;
plot(snr,RMSE_elev_E_c_GMNS,'m-s');
legend('D-SHESPRIT-GMNS');
xlabel('SNR (dB)');
ylabel('RMSE of elevation (degree)');

figure;
plot(snr,RMSE_azim_E_c_GMNS,'m-s');
legend('D-SHESPRIT-GMNS');
xlabel('SNR (dB)');
ylabel('RMSE of azimuth (degree)');