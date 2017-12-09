%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Demo-Script Structured Least Squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script demonstrates the basic version of Structured Least Squares
% (one iteration, no regularization) for a Uniform Linear Array in
% conjunction with 1-D standard ESPRIT. It is easy to adapt for (R-D)
% Unitary ESPRIT and nonuniform array geometries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Florian Roemer, DVT, TU Ilmenau, May 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%% ULA setup
% number of array elements
M = 15; 
% array manifold
A_mu = @(mu) exp(1i*(0:M-1)'*mu);
% selection matrices for ESPRIT
J1 = speye(M-1,M);
J2 = rot90(J1,2);

%% Scenario
% source positions
mu0 = [0.3,1.5, -1.2];
% number of snapshots
T = 3;
% signal to noise ratio [dB]
SNR = 20;

%% Data model
% number of sources
K = length(mu0);
% amplitudes (here: complex Gaussian)
S = (randn(K,T)+1i*randn(K,T))/sqrt(2);
% noise (here: AWGN ~ZMCSCG)
W = (randn(M,T)+1i*randn(M,T))/sqrt(2)*10^(-SNR/20);
% data
X = A_mu(mu0)*S + W;

%% Subspace estimation
% model order (here: assumed perfectly known)
Khat = K;
% singular value decomposition
[U,~] = svd(X);
% subspace: dominant LSV
Us = U(:,1:Khat);

%% LS Standard ESPRIT
% solve shift invariance equation (LS)
Psi_LS = (J1*Us)\(J2*Us);    
% estimate angles
muhat_LS = angle(eig(Psi_LS));

%% SLS Standard ESPRIT
% residual (from LS solution)
R_LS = J1*Us*Psi_LS - J2*Us;
% helper matrices (F = [IJ1Us,PsiJ12])
IJ1Us = kron(speye(Khat),J1*Us);
PsiJ12 = kron(Psi_LS.',J1)-kron(speye(Khat),J2);   
% sls update (vectorized, upper part of pinv(F)*r_LS)
upd_sls = -IJ1Us' * (( IJ1Us*IJ1Us' + PsiJ12*PsiJ12' ) \ R_LS(:));
% apply update to LS solution
Psi_SLS = Psi_LS + reshape(upd_sls,[Khat,Khat]);
% estimate angles
muhat_SLS = angle(eig(Psi_SLS));

%% Evaluate mean square estimation error
% here: quick'n'dirty via sort (better: use association algorithm, e.g.,
% Hungarian algorithm)
MSE_LS = sum((sort(mu0(:),'descend')-sort(muhat_LS,'descend')).^2);
MSE_SLS = sum((sort(mu0(:),'descend')-sort(muhat_SLS,'descend')).^2);
% display results
fprintf('MSE  LS = %g\nMSE SLS = %g\n',MSE_LS,MSE_SLS);