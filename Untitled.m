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
for ii=1:length(snr)
    refp=10.^(snr(ii)/20);
    S=refp.*ss;
    [q,T]=size(S);
     Signal=A*S;
for p=1:monte_num
  
         noise=randn(M,snapshot)+j*randn(M,snapshot);

         x=Signal+noise;
         C=pinv(y_omega);%tubu
        C2=pinv(Bn);
         x_y=C2*C*x; 
         tic
         Rnm=cov(x_y');
        [U,V]=eig(Rnm);        
         Sspace=U(:,(L1-I+1):L1);
      toc
       disp(['evd运行时间',num2str(toc)]);
          tic
        for i = 1:K
       X_re = x_y((i-1)*l+1:i * l,:);
       R_re = 1/T * (X_re * X_re');
       [Vec,Val] = eig(R_re);
       [~, Pos1] = sort(diag(Val),'descend');
       W_re = Vec(:,Pos1(1:q));
       S_1 = pinv(W_re) * X_re;
       if i ==1
           W = W_re.';
           S1 = S_1;
           Q=(trace(R_re)-trace(W_re'*R_re*W_re))/(l-q);
       else
           T_est = (S_1*S1'/T)*(S1*S1'/T-eye(q)*Q)^-1;
           W = [W,(W_re * T_est).'];
       end
   end
   W = W.';
   Sspace_GMNS=W;
toc
disp(['gmns运行时间',num2str(toc)]);
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

      elev_estimated_c(:,p)=acos(diag(Dm_c))*180/pi;
      elev_estimated_c=sort(elev_estimated_c);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT GMNS %%%%%%%%%%%%%%%%%%%%%%%% 
      SS1_c_GMNS=pinv(S1_c_GMNS'*S1_c_GMNS);
      FAI_c_GMNS=SS1_c_GMNS*S1_c_GMNS'*(La_c_GMNS*S0_c_GMNS+La2_c_GMNS*S2_c_GMNS);
      [Vm_c_GMNS,Dm_c_GMNS]=eig(FAI_c_GMNS);

      elev_estimated_c_GMNS(:,p)=acos(diag(Dm_c_GMNS))*180/pi;
      elev_estimated_c_GMNS=sort(elev_estimated_c_GMNS);
      %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%%  
            SS_a11_c=pinv(S_a0_c'*S_a0_c);
      FAI_a1_c=SS_a11_c*S_a0_c'*(-Mu1_c*S_a1_c+Mu2_c*S_a2_c);
      [Vm_c1,Dm_c1]=eig(FAI_a1_c);

      azim_estimated_c(:,p)=angle(diag(Dm_c1)./(sin(acos(diag(Dm_c)))))*180/pi;  
      azim_estimated_c=sort(azim_estimated_c);    
       %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT GMNS %%%%%%%%%%%%%%%%%%%%%%%%%%      
      SS_a11_c_GMNS=pinv(S_a0_c_GMNS'*S_a0_c_GMNS);
      FAI_a1_c_GMNS=SS_a11_c_GMNS*S_a0_c_GMNS'*(-Mu1_c_GMNS*S_a1_c_GMNS+Mu2_c_GMNS*S_a2_c_GMNS);
      [Vm_c1,Dm_c1_GMNS]=eig(FAI_a1_c_GMNS);

      azim_estimated_c_GMNS(:,p)=angle(diag(Dm_c1_GMNS)./(sin(acos(diag(Dm_c_GMNS)))))*180/pi;  
      azim_estimated_c_GMNS=sort(azim_estimated_c_GMNS);      
end
     RMSE_azim_esprit_c_GMNS(:,ii)=sqrt(sum((((azim_estimated_c_GMNS-azimuth*ones(1,p))).^2)')/p);
     RMSE_elev_estimated_c_GMNS(:,ii)=sqrt(sum((abs(((elev_estimated_c_GMNS-elevation*ones(1,p)))).^2)')/p);
     
     RMSE_azim_esprit_c(:,ii)=sqrt(sum((((azim_estimated_c-azimuth*ones(1,p))).^2)')/p);
     RMSE_elev_estimated_c(:,ii)=sqrt(sum((abs(((elev_estimated_c-elevation*ones(1,p)))).^2)')/p);
end

RMSE_azim_E_c_GMNS=sum(RMSE_azim_esprit_c_GMNS)/2;
RMSE_elev_E_c_GMNS=sum(RMSE_elev_estimated_c_GMNS)/2;

RMSE_azim_E_c=sum(RMSE_azim_esprit_c)/2;
RMSE_elev_E_c=sum(RMSE_elev_estimated_c)/2;

figure;
plot(snr,RMSE_elev_E_c_GMNS,'m-s',snr,RMSE_elev_E_c,'b-^');
legend('D-SHESPRIT-GMNS','D-SHESPRIT');
xlabel('SNR (dB)');
ylabel('RMSE of elevation (degree)');

figure;
plot(snr,RMSE_azim_E_c_GMNS,'m-s',snr,RMSE_azim_E_c,'b-^');
legend('D-SHESPRIT-GMNS','D-SHESPRIT');
xlabel('SNR (dB)');
ylabel('RMSE of azimuth (degree)');
