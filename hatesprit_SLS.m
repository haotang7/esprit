function DOA=hatesprit_SLS(Us,Fai,noDOA)
%***********************************************************************
%author:Shiyunmei
%time:2013.05.29 22:14
%input:
%Xx: the received signal
%noDOA: number of DOA
%output:
%DOA: the estimated angle of signal
%***********************************************************************
% [M,N]=size(Xx);
% M=M/m_time;
aefa=3;
epsion=1e-3;
% Ka=0.1;
% Ka=sqrt(9/(aefa*25));
Ka=sqrt(noDOA/(aefa*25*9));

S1_r1=Us([3,8,7,6,15,14,13,12,11],:);
S0_r1=Us([1,4,3,2,9,8,7,6,5],:);
S2_r1=Us([7,14,13,12,23,22,21,20,19],:);
        
t1_r1=zeros(9,25);t0_r1=zeros(9,25);t2_r1=zeros(9,25);
t1_r1(1,3)=1;t1_r1(2,8)=1;t1_r1(3,7)=1;t1_r1(4,6)=1;t1_r1(5,15)=1;t1_r1(6,14)=1;t1_r1(7,13)=1;t1_r1(8,12)=1;t1_r1(9,11)=1;
t0_r1(1,1)=1;t0_r1(2,4)=1;t0_r1(3,3)=1;t0_r1(4,2)=1;t0_r1(5,9)=1;t0_r1(6,8)=1;t0_r1(7,7)=1;t0_r1(8,6)=1;t0_r1(9,5)=1;
t2_r1(1,7)=1;t2_r1(1,14)=1;t2_r1(1,13)=1;t2_r1(1,12)=1;t2_r1(1,23)=1;t2_r1(1,22)=1;t2_r1(1,21)=1;t2_r1(1,20)=1;t2_r1(1,19)=1;

lamda1_r1=[lamda(1,0),lamda(2,-1),lamda(2,0),lamda(2,1),lamda(3,-2),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(3,2)];
lamda2_r1=[lamda(2,0),lamda(3,-1),lamda(3,0),lamda(3,1),lamda(4,-2),lamda(4,-1),lamda(4,0),lamda(4,1),lamda(4,2)];
La11=diag(lamda1_r1);
La21=diag(lamda2_r1);

tt3=La11*t0_r1+La21*t2_r1;
% Rx=Xx*Xx'/N;
% [U,S] = svd(Rx);
% Us = U(:,1:noDOA);
% J1=kron(eye(m_time),[eye(M-1),zeros(M-1,1)]);
% J2=kron(eye(m_time),[zeros(M-1,1),eye(M-1)]);
% Us1=J1*Us;
% Us2=J2*Us;
% Psi=Us1\Us2;
% U=Us;%the initial basis of the estimated signal subspace for k=1


delta_Usk=zeros(size(Us));
Id=eye(noDOA);
Z=[kron(Id,S1_r1) kron(Fai.',t1_r1)-kron(Id,tt3);zeros(25*noDOA,noDOA*noDOA),Ka*eye(25*noDOA)];
R=S1_r1*Fai-La11*S0_r1-La21*S2_r1;
b_re=[vec(R);Ka*vec(delta_Usk)];
delta_PsiU=-pinv(Z'*Z)*Z'*b_re;
delta_Psi=delta_PsiU(1:noDOA*noDOA);
delta_U=delta_PsiU(noDOA*noDOA+1:end);
% while norm(delta_Psi)^2>epsion && norm(delta_U)^2>epsion
while max(norm(delta_Psi)^2,norm(delta_U)^2)>epsion
    delta_Usk=delta_Usk+reshape(delta_U,25,noDOA);
    Us=Us+reshape(delta_U,25,noDOA);
    Fai=Fai+reshape(delta_Psi,noDOA,noDOA);
    Z=[kron(Id,t1_r1*Us) kron(Fai.',t1_r1)-kron(Id,tt3);zeros(25*noDOA,noDOA*noDOA),Ka*eye(25*noDOA)];
    R=t1_r1*Us*Fai-La11*t0_r1*Us-La21*t2_r1*Us;
    b_re=[vec(R);Ka*vec(delta_Usk)];
    delta_PsiU=-pinv(Z'*Z)*Z'*b_re;
    delta_Psi=delta_PsiU(1:noDOA*noDOA);
    delta_U=delta_PsiU(noDOA*noDOA+1:end);
end
[~,D] = eig(Fai);
DOA = acos(diag(D))*180/pi;
end