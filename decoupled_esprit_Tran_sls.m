% close all;
clear variables;
tic
j=sqrt(-1);
% ======= elevation and azimuth in degree ======== %
azimuth=[42;149];      %方位角[60;128] 
elevation=[55;78];   %俯仰角 [49;70]  
degrad=pi/180;%角度-弧度换算
azimuth1=azimuth*degrad;
elevation1=elevation*degrad;
% ===================signal and spherical array============================= %
snapshot=200;%快拍数
N_x=(1:snapshot);%采样点数；
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
L1=(N+1).^2;%Y矩阵的列长度
%%%%%%%%---------------------求出A矩阵--------------------------------%%%%%%%
y_omega_r=zeros(M,L1);%生成零矩阵Y(omega)
y_fai_r=zeros(1,L1);%生成零矩阵Y(fai)
y_fai_r_me=zeros(1,L1);%生成零矩阵Y(fai)
y_omega=zeros(M,L1);%生成零矩阵Y(omega)
y_fai=zeros(1,L1);%生成零矩阵Y(fai)

for n=0:N
    for m=-n:n 
        for  k1=1:M
             y_omega_r(k1,n.^2+n+m+1)=spherical_har_real(m,n,theta_y(k1),fai_y(k1));     % 生成Y(omega)矩阵
        end
        for i=1:I                
            y_fai_r(i,n.^2+n+m+1)=spherical_har_real(m,n,elevation1(i),azimuth1(i));     % 生成Y(fai)矩阵
        end

    end
end 

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


Y_f_r=y_fai_r';
Y_fai_r=Bn*y_fai_r';%求Y（fai）的复共轭
A_r=y_omega_r* Y_fai_r;


Y_f=y_fai';
Y_fai=Bn*y_fai';%求Y（fai）的复共轭
A=y_omega* Y_fai;


ss=randn(I,snapshot)+j*randn(I,snapshot);

monte_num=222;
snr=0:5:20;
Tran=Trans3(25);

for ii=1:length(snr)
    refp=10.^(snr(ii)/20);
    S=refp.*ss;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Rs=cov(S');
     Signal=A*S;
     p1=0;
     azi_estimation=[];
     ele_estimation=[];
      Doa3st_3=[];
      Doa3st_3_r=[];
      azi_estimation_m=[];
      azi_estimation_m1=[];
            azim_m1=[];
%       crb=[];
    p2 = 0;
    p3 = 0;
    p4=0;
for p=1:monte_num
  
         noise=randn(M,snapshot)+j*randn(M,snapshot);

         x=Signal+noise;
        C=pinv(y_omega);%tubu
        C2=pinv(Bn);
         x_y=C2*C*x;       
         Rnm=cov(x_y');
        [U,V]=eig(Rnm);

      %%%%%%%%%%%%% 只乘Q %%%%%%%%%%%%%%%%%%%%  
         x_y_r2=Tran*x_y;
         Rnm_r2=cov(x_y_r2');
         [U0_r2,V0_r2]=eig(Rnm_r2);

      %%%%%%%%%%%%%%%%实虚分开%%%%%%%%%%%%%
         x_y_r1=[real(x_y_r2),imag(x_y_r2)];
         Rnm_r1=cov(x_y_r1');
         [U0_r1,V0_r1]=eig(Rnm_r1);
         
 %%%%%%%%%%%%%%%%%%%%%%%% real spherical harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Sspace_r1=U0_r1(:,(L1-I+1):L1);%%%% RV-SHESPRIT
        Sspace_r2=U0_r2(:,(L1-I+1):L1);%%% U-SHESPRIT
        
        S1_r1=Sspace_r1([3,8,7,6,15,14,13,12,11],:);
        S0_r1=Sspace_r1([1,4,3,2,9,8,7,6,5],:);
        S2_r1=Sspace_r1([7,14,13,12,23,22,21,20,19],:);
           lamda1_r1=[lamda(1,0),lamda(2,-1),lamda(2,0),lamda(2,1),lamda(3,-2),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(3,2)];
           lamda2_r1=[lamda(2,0),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(4,-2),lamda(4,-1),lamda(4,0),lamda(4,1),lamda(4,2)];
           La11=diag(lamda1_r1);
           La21=diag(lamda2_r1);

           
        S_a01=Sspace_r1([2,6,5,12,11,10],:);
        S_a11=Sspace_r1([4,8,9,14,15,16],:);
        S_a21=Sspace_r1([1,3,2,7,6,5],:);
        S_a31=Sspace_r1([7,13,12,21,20,19],:);
        
                
           mu1_r=[mu(1,1),mu(2,1),mu(2,2),mu(3,1),mu(3,2),mu(3,3)];
           mu2_r=[mu(2,0),mu(3,0),mu(3,-1),mu(4,0),mu(4,-1),mu(4,-2)];
           mu3_r=[sqrt(2),sqrt(2),1,sqrt(2),1,1];
           Mu1_r=diag(mu1_r);
           Mu2_r=diag(mu2_r);
           Mu3_r=diag(mu3_r);
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           S1_r2=Sspace_r2([3,8,7,6,15,14,13,12,11],:);
           S0_r2=Sspace_r2([1,4,3,2,9,8,7,6,5],:);
           S2_r2=Sspace_r2([7,14,13,12,23,22,21,20,19],:);

           
           S_a02=Sspace_r2([2,6,5,12,11,10],:);
           S_a12=Sspace_r2([4,8,9,14,15,16],:);
           S_a22=Sspace_r2([1,3,2,7,6,5],:);
           S_a32=Sspace_r2([7,13,12,21,20,19],:);
        

        %%%%%%%%%%%%%%%%%%% complex spherical harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Sspace=U(:,(L1-I+1):L1);
         
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

            t1_r1=zeros(9,L1);t0_r1=zeros(9,L1);t2_r1=zeros(9,L1);   
            
            t1_r1(1,3)=1;t1_r1(2,6)=1;t1_r1(3,7)=1;t1_r1(4,8)=1;t1_r1(5,11)=1;t1_r1(6,12)=1;t1_r1(7,13)=1;t1_r1(8,14)=1;t1_r1(9,15)=1;
           t0_r1(1,1)=1;t0_r1(2,2)=1;t0_r1(3,3)=1;t0_r1(4,4)=1;t0_r1(5,5)=1;t0_r1(6,6)=1;t0_r1(7,7)=1;t0_r1(8,8)=1;t0_r1(9,9)=1;
           t2_r1(1,7)=1;t2_r1(1,12)=1;t2_r1(1,13)=1;t2_r1(1,14)=1;t2_r1(1,19)=1;t2_r1(1,20)=1;t2_r1(1,21)=1;t2_r1(1,22)=1;t2_r1(1,23)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%
 S_1=[];S0=[];S1=[];tn=[];r1=zeros(N^2,N^2);r2=zeros(N^2,N^2);
        for n=1:N
            % S_1:抽取前N-1个元素
            n_temp1=n^2+1;
            s_temp1=Sspace(n_temp1:(n_temp1+2*n-2),:);
            S_1=[S_1;s_temp1];
             
        % S0:抽取中间N-1个元素
            n_temp2=n^2+2;
            s_temp2=Sspace(n_temp2:(n_temp2+2*n-2),:);
            S0=[S0;s_temp2];
             
        % S1:抽取后N-1个元素
            n_temp3=n^2+3;
            s_temp3=Sspace(n_temp3:(n_temp3+2*n-2),:);
            S1=[S1;s_temp3];
             
        %------------------- 构建抽取矩阵的系数矩阵 --------------------%
            if n==1
                t_temp=0;
            else
                t_temp=(-(n-1):1:(n-1));
            end
            tn=[tn,t_temp];
             
            for m=-(n-1):(n-1)
                lamda1=sqrt((n-m)*(n+m+1));
                lamda2=sqrt((n+m)*(n-m+1));
                r1((n-1)^2+(n-1)+m+1,(n-1)^2+(n-1)+m+1)=lamda1;
                r2((n-1)^2+(n-1)+m+1,(n-1)^2+(n-1)+m+1)=lamda2;
            end  
        end

        T=diag(tn);

    % ============================ 计算声源的距离信息 ================================= %
        E=[r2*S_1,r1*S1];
        FAI=-2*inv((E'*E))*E'*T*S0;        % FAI是一2I*I矩阵
        FQ=FAI(1:I,:);                     % 抽取前I*I矩阵
        [Vm,Dm]=eig(FQ);                   % 特征分解，特征值对应声源信息
        
    % compute angles
        azim_estimated(:,p)=-angle(diag(Dm))*180/pi;
        elev_estimated(:,p)=atan(abs(diag(Dm)))*180/pi;
        
        azim_estimated=sort(azim_estimated);
        elev_estimated=sort(elev_estimated);
 
 %%%%%%%%%%%%%%%%%%%%%%%% D-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%

      SS1_c=pinv(S1_c'*S1_c);
      FAI_c=SS1_c*S1_c'*(La_c*S0_c+La2_c*S2_c);
      [Vm_c,Dm_c]=eig(FAI_c);

      elev_estimated_c(:,p)=acos(diag(Dm_c))*180/pi;
      elev_estimated_c=sort(elev_estimated_c);
      
  %%%%%%%%%%%%%%%%%%%%% SLS-D-SHESPRIT &&&&&&&&&&&&&&&&&&&&&&&&&&
  
        RR=S1_c*FAI_c-La_c*S0_c-La2_c*S2_c;
        Z11=kron(speye(I),S1_c);
        tt3=La_c*t0_r1+La2_c*t2_r1;
        Z12=kron(FAI_c.',t1_r1)-kron(speye(I),tt3);
        upd_sls=-Z11'*((Z11*Z11'+Z12*Z12')\RR(:));
        FAI_sls=FAI_c+reshape(upd_sls,[I,I]);
        
        [Vm_sls,Dm_sls]=eig(FAI_sls);

      elev_estimated_sls(:,p)=acos(diag(Dm_sls))*180/pi;
      elev_estimated_sls=sort(elev_estimated_sls);
      
      
      
    %%%%%%%%%%%%%%%%%%%%%%%% New-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%%%

      SS_a11_c=pinv(S_a0_c'*S_a0_c);
      FAI_a1_c=SS_a11_c*S_a0_c'*(-Mu1_c*S_a1_c+Mu2_c*S_a2_c);
      [Vm_c1,Dm_c1]=eig(FAI_a1_c);

      azim_estimated_c(:,p)=angle(diag(Dm_c1))*180/pi;
      elev_estimated_c2(:,p)=asin(abs(diag(Dm_c1)))*180/pi;
      
      elev_estimated_c2=sort(elev_estimated_c2);
      azim_estimated_c=sort(azim_estimated_c);

     %%%%%%%%%%%%%%%%%%%%%%%% RV-SHESPRIT %%%%%%%%%%%%%%%%%%%%%%%%
    
      SS1_r1=pinv(S1_r1'*S1_r1);
      FAI_r1=SS1_r1*S1_r1'*(La11*S0_r1+La21*S2_r1);
      [Vm_r1,Dm_r1]=eig(FAI_r1);
      TZ1=inv(Vm_r1);

      elev_estimated_r1(:,p)=acos(diag(Dm_r1))*180/pi;
      elev_estimated_r1=sort(elev_estimated_r1);
      
       %%%%%%%%%%%%%%%%%%%%% SLS-RV-SHESPRIT &&&&&&&&&&&&&&&&&&&&&&&&&&
  
        RR_r1=S1_r1*FAI_r1-La11*S0_r1-La21*S2_r1;
        Z11_r1=kron(speye(I),S1_r1);
        tt3_r1=La11*t0_r1+La21*t2_r1;
        Z12_r1=kron(FAI_r1.',t1_r1)-kron(speye(I),tt3_r1);
        upd_sls_r1=-Z11_r1'*((Z11_r1*Z11_r1'+Z12_r1*Z12_r1')\RR_r1(:));
        FAI_sls_r1=FAI_r1+reshape(upd_sls_r1,[I,I]);
        
        [Vm_sls_r1,Dm_sls_r1]=eig(FAI_sls_r1);

      elev_estimated_sls_r1(:,p)=acos(diag(Dm_sls_r1))*180/pi;
      elev_estimated_sls_r1=sort(elev_estimated_sls_r1);
      
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      T_n1=blkdiag(pinv(TZ1),pinv(TZ1));
      Z1=TZ1;

      SL1=[S_a01,S_a11];
      SS_a11=pinv(SL1'*SL1);
      FAI_a1=-SS_a11*SL1'*Mu3_r*(Mu1_r*S_a21-Mu2_r*S_a31);
      TTT_a1=pinv(T_n1)*FAI_a1*pinv(Z1);
             
      FAI_a111=TTT_a1(1:2,:);
      FAI_a221=TTT_a1(3:4,:);

      FAI11=FAI_a221(1,1)/FAI_a111(1,1);
      FAI21=FAI_a221(2,2)/FAI_a111(2,2);
      FAI31=[FAI11,FAI21];
      FAI41=diag(FAI31);

      azim_estimated_r1(:,p)=atan(diag(FAI41))*180/pi;
     if azim_estimated_r1(1,p)<0
          azim_estimated_r1(1,p)=180-abs(azim_estimated_r1(1,p));
      end
      if azim_estimated_r1(2,p)<0
          azim_estimated_r1(2,p)=180-abs(azim_estimated_r1(2,p));
      end

      azim_estimated_r1=sort(azim_estimated_r1);
      

     
end
   
     RMSE_azim_esprit(:,ii)=sqrt(sum((((azim_estimated-azimuth*ones(1,p))).^2)')/p);
     RMSE_azim_esprit_c(:,ii)=sqrt(sum((((azim_estimated_c-azimuth*ones(1,p))).^2)')/p);
%      RMSE_azim_esprit_c2(:,ii)=sqrt(sum(((azim_estimated_c2-azimuth*ones(1,p)).^2)')/p);
     RMSE_azim_estimated_r1(:,ii)=sqrt(sum((((azim_estimated_r1-azimuth*ones(1,p))).^2)')/p);
%      RMSE_azim_estimated_r2(:,ii)=sqrt(sum((abs((azim_estimated_r2-azimuth*ones(1,p))).^2)')/p);
     
     RMSE_elev_esprit(:,ii)=sqrt(sum((((elev_estimated-elevation*ones(1,p))).^2)')/p);
     RMSE_elev_esprit_c(:,ii)=sqrt(sum((abs(((elev_estimated_c-elevation*ones(1,p)))).^2)')/p);
      RMSE_elev_esprit_c2(:,ii)=sqrt(sum((((elev_estimated_c2-elevation*ones(1,p))).^2)')/p);
     RMSE_elev_estimated_r1(:,ii)=sqrt(sum((((elev_estimated_r1-elevation*ones(1,p))).^2)')/p);
     RMSE_elev_SLS_estimated_c(:,ii)=sqrt(sum((((elev_estimated_sls-elevation*ones(1,p))).^2)')/p);
     RMSE_elev_SLS_estimated_r1(:,ii)=sqrt(sum((((elev_estimated_sls_r1-elevation*ones(1,p))).^2)')/p);
%      RMSE_elev_estimated_r2(:,ii)=sqrt(sum((abs((elev_estimated_r2-elevation*ones(1,p))).^2)')/p);

     %%%%%%%%%%%%%%%%%%%%%%%%
   crb1(:,ii)=ball_crb(M,N,elevation1,azimuth1,y_omega,Bn,Rs,snapshot,1);      
   crb_ele1(:,ii)=crb1(1,ii);
   crb_ele2(:,ii)=crb1(2,ii);
   crb_ele(:,ii)=(crb_ele1(:,ii)+crb_ele2(:,ii))/2;
   
   crb_azi1(:,ii)=crb1(3,ii);
   crb_azi2(:,ii)=crb1(4,ii);
   crb_azi(:,ii)=(crb_azi1(:,ii)+crb_azi2(:,ii))/2;

end

RMSE_azim_E=sum(RMSE_azim_esprit)/2;
RMSE_azim_E_c=sum(RMSE_azim_esprit_c)/2;
% RMSE_azim_E_c2=sum(RMSE_azim_esprit_c2)/2;
RMSE_azim_E_r1=sum(RMSE_azim_estimated_r1)/2;
% RMSE_azim_E_r2=sum(RMSE_azim_estimated_r2)/2;

RMSE_elev_E=sum(RMSE_elev_esprit)/2;
RMSE_elev_E_c=sum(RMSE_elev_esprit_c)/2;
RMSE_elev_E_c2=sum(RMSE_elev_esprit_c2)/2;
RMSE_elev_E_r1=sum(RMSE_elev_estimated_r1)/2;
RMSE_elev_SLS_E_c=sum(RMSE_elev_SLS_estimated_c)/2;
RMSE_elev_SLS_E_r1=sum(RMSE_elev_SLS_estimated_r1)/2;
% RMSE_elev_E_r2=sum(RMSE_elev_estimated_r2)/2;


% 
figure;
% set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',2);
% set(gcf,'Position',[347,162,330,300]);
semilogy(snr,RMSE_elev_E,'r-o',snr,RMSE_elev_E_c,'c-s',snr,RMSE_elev_SLS_E_c,'m-h',snr,RMSE_elev_E_r1,'b-^',snr,RMSE_elev_SLS_E_r1,'g-h',snr,crb_ele,'k--');
legend('SHESPRIT','D-SHESPRIT','SLS-D-SHESPRIT','RV-SHESPRIT','SLS-RV-SHESPRIT','CRB');
xlabel('SNR (dB)');
ylabel('RMSE of elevation (degree)');

figure;
% set(gca,'FontName','Times New Roman','FontSize',9,'LineWidth',2);
% set(gcf,'Position',[347,162,330,300]);
plot(snr,RMSE_azim_E,'r-o',snr,RMSE_azim_E_c,'c-s',snr,RMSE_azim_E_r1,'b-^',snr,crb_azi,'k--');
legend('SHESPRIT','D-SHESPRIT','RV-SHESPRIT','CRB');
xlabel('SNR (dB)');
ylabel('RMSE of azimuth (degree)');



  