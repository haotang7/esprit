clear variables;%清除变量
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
N_x=(1:snapshot);%采样点数
I=size( azimuth, 1 );% number of  signals
load mic_position;
theta_y=mic_position(:,1)*degrad;
fai_y=mic_position(:,2)*degrad;
M=32;
N=4;
f=[150 200];    %信号调制频率
fs=4500;            %采样频率(满足采样定理)
c=343;
r=0.1;
kr=[4;4];
ka=kr;
L1=(N+1).^2;%Y矩阵的列长度p33
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
Y_f=y_fai';%'zhuanzhi
Y_fai=Bn*y_fai';%求Y（fai）的复共轭
A=y_omega* Y_fai;
ss=randn(I,snapshot)+j*randn(I,snapshot);
monte_num=100;
snr=0:5:20;     %%对信噪比离散化取值
Tran=Trans3(25);
for ii=1:length(snr)
    refp=10.^(snr(ii)/20);
    S=refp.*ss;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Signal=A*S;
for p=1:monte_num
        noise=randn(M,snapshot)+j*randn(M,snapshot);   
        x=Signal+noise;
        C=pinv(y_omega);%tubu
        C2=pinv(Bn);
        x_y=C2*C*x;
        Rnm=cov(x_y');%求方差
        [U,V]=eig(Rnm);   %[V,D]=eig(A)：求矩阵A的全部特征值，构成对角阵D，并求A的特征向量构成V的列向量    
         Sspace=U(:,(L1-I+1):L1);
%%%%%%%%%%%%%%%%%%%%%%%%% complex spherical harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S1_c=Sspace([3,6,7,8,11,12,13,14,15],:);
        S0_c=Sspace([1,2,3,4,5,6,7,8,9],:);
        S2_c=Sspace([7,12,13,14,19,20,21,22,23],:);
           lamda1_c=[lamda(1,0),lamda(2,-1),lamda(2,0),lamda(2,1),lamda(3,-2),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(3,2)];
           lamda2_c=[lamda(2,0),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(4,-2),lamda(4,-1),lamda(4,0),lamda(4,1),lamda(4,2)];
           La_c=diag(lamda1_c);
           La2_c=diag(lamda2_c);
           t1_c=zeros(9,L1);t0_c=zeros(9,L1);t2_c=zeros(9,L1);        
           t1_c(1,3)=1;t1_c(2,6)=1;t1_c(3,7)=1;t1_c(4,8)=1;t1_c(5,11)=1;t1_c(6,12)=1;t1_c(7,13)=1;t1_c(8,14)=1;t1_c(9,15)=1;
           t0_c(1,1)=1;t0_c(2,2)=1;t0_c(3,3)=1;t0_c(4,4)=1;t0_c(5,5)=1;t0_c(6,6)=1;t0_c(7,7)=1;t0_c(8,8)=1;t0_c(9,9)=1;
           t2_c(1,7)=1;t2_c(2,12)=1;t2_c(3,13)=1;t2_c(4,14)=1;t2_c(5,19)=1;t2_c(6,20)=1;t2_c(7,21)=1;t2_c(8,22)=1;t2_c(9,23)=1;   
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
           S_a0_c=Sspace([4,8,9,12,13,14,15,16],:);
           S_a1_c=Sspace([1,3,4,5,6,7,8,9],:);
           S_a2_c=Sspace([7,13,14,19,20,21,22,23],:);
           mu1_c=[mu(1,1),mu(2,1),mu(2,2),mu(3,-1),mu(3,0),mu(3,1),mu(3,2),mu(3,3)];
           mu2_c=[mu(2,0),mu(3,0),mu(3,-1),mu(4,2),mu(4,1),mu(4,0),mu(4,-1),mu(4,-2)];          
           Mu1_c=diag(mu1_c);
           Mu2_c=diag(mu2_c);
           tt1_c=zeros(8,L1);tt0_c=zeros(8,L1);tt2_c=zeros(8,L1);   
           tt0_c(1,4)=1;tt0_c(2,8)=1;tt0_c(3,9)=1;tt0_c(4,12)=1;tt0_c(5,13)=1;tt0_c(6,14)=1;tt0_c(7,15)=1;tt0_c(8,16)=1;
           tt1_c(1,1)=1;tt1_c(2,3)=1;tt1_c(3,4)=1;tt1_c(4,5)=1;tt1_c(5,6)=1;tt1_c(6,7)=1;tt1_c(7,8)=1;tt1_c(8,9)=1;
           tt2_c(1,7)=1;tt2_c(2,13)=1;tt2_c(3,14)=1;tt2_c(4,19)=1;tt2_c(5,20)=1;tt2_c(6,21)=1;tt2_c(7,22)=1;tt2_c(8,23)=1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%
      SS1_c=pinv(S1_c'*S1_c);
      FAI_c=SS1_c*S1_c'*(La_c*S0_c+La2_c*S2_c);
      [Vm_c,Dm_c]=eig(FAI_c);
      elev_estimated_c(:,p)=acos(diag(Dm_c))*180/pi;
      elev_estimated_c=sort(elev_estimated_c);
      %%%%%%%%%%%%%%%%%%%% Azimuth D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%% 
      SS_a11_c=pinv(S_a0_c'*S_a0_c);
      FAI_a1_c=SS_a11_c*S_a0_c'*(-Mu1_c*S_a1_c+Mu2_c*S_a2_c);
      [Vm_c1,Dm_c1]=eig(FAI_a1_c);
      azim_estimated_c(:,p)=angle(diag(Dm_c1)./(sin(acos(diag(Dm_c)))))*180/pi;  
      azim_estimated_c=sort(azim_estimated_c);
end
     RMSE_azim_esprit_c(:,ii)=sqrt(sum((((azim_estimated_c-azimuth*ones(1,p))).^2)')/p);
     RMSE_elev_estimated_c(:,ii)=sqrt(sum((abs(((elev_estimated_c-elevation*ones(1,p)))).^2)')/p);
end
RMSE_azim_E_c=sum(RMSE_azim_esprit_c)/2;
RMSE_elev_E_c=sum(RMSE_elev_estimated_c)/2;
figure;
plot(snr,RMSE_elev_E_c,'m-s',snr,RMSE_azim_E_c);
legend('elevation','azimuth');
xlabel('SNR (dB)');
ylabel('RMSE of(degree)');




