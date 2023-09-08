% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%
clear variables
warning('on','all')
rnd_seed = 10; %randi(1000);

% This example shows how the lib can be used to simulate the SER in
% an FSC channel
% Notes: 1) The theoretical SER are still in progress
%        2) The proper channel equalization approach is a open topic to
%        work with

%% Simulation Parameters

% Estimation techniques
white_LMMSE_Zm = true;
detection_flag = true;
HPC_flag = false;


if HPC_flag
    MCS = input('Please choose the MCS \n (1) QPSK 1/3 \n (2) 16-QAM 1/3 \n (3) 16-QAM 2/3 \n (4) 64-QAM 3/4. \n');
    fd = input('Doppler shift? ');
    its = input('Number of simulation trials? ');
    rnd_seed = input('seed number: ');   

    parpool(16)    

    CurrentFolder = pwd
%    addpath(genpath(CurrentFolder))
    
%   mex ./Lib/cml/source/ConvEncode.c
%   mex ./Lib/cml/source/Deinterleave.c
%   mex ./Lib/cml/source/Demod2D.c
%   mex ./Lib/cml/source/Depuncture.c
%   mex ./Lib/cml/source/Interleave.c
%   mex ./Lib/cml/source/Modulate.c
%   mex ./Lib/cml/source/Puncture.c
%   mex ./Lib/cml/source/Somap.c
%   mex ./Lib/cml/source/SisoDecode.c


else
    MCS = 3;
    fd = 1500;%0.1;
    its = 10;
    rnd_seed = 200;
    
    % TEST SNR 
    % (in case of normal simulation must be commented out)
    conf.test_snr = 10;
end


max_ce_input_snr = 10; % Maximumum SNR for the input of ZAK based channel estimators
max_ce_input_snr_high = 30; 
%snr = [-5:5:35 60];


% snr = [20:5:35 60];
p=get_defaultGFDM('OFDM');
p.mode = '2blocks';
p.pulse='rc_fd';
p.a=0.3;%roll-off factor for the RRC and RC filters.
p.K=3*48;%Number of useful subcarriers.
%p.L=2;%Number of overlapping subcarrier for frequency domain GFDM implementation..
p.Mp = 2; % Number of pilot subsymbols
p.Md = 5; p.M= p.Md + p.Mp; % Number of data subsymbols %Number of subsymbols.
%p.Md = 7; p.M= p.Md; % Number of data subsymbols %Number of subsymbols.

%Koff = 11;  % Total number of zero subcarriers
Koff = 3;

k_start = ceil(Koff/2);
k_end = floor(Koff/2);
p.Kset = (k_start+1:(p.K-k_end) )-1;
p.Kon=length(p.Kset);%Number of actives subcarriers.
p.k_start = k_start;
p.k_end = k_end;
if MCS == 1
   % snr = [-5 -2 0 3 5 10] %[-5:5:35]; % [0 10 20 30 50 80 100 200]; %0:5:35; % 0:10:50;
    snr = [-5 -2 0 3 5 7 9 11]
  %  high_snr = 5;
    %snr = [-4.5 -3.5 -3 -2.75 -1 5 10 20]
    p.mu = 2;
    p.coderate = 1/3; % 1/3; Puncturing pattern is defined for 1/3 & 3/4 only
elseif MCS == 2
    snr = [-5 0 5 10 12 14 15]
   % high_snr = 12;
    p.mu = 4;
    p.coderate = 1/3;
elseif MCS == 3 
    snr = [5 10 15 20 25 30 35 40]
   % high_snr = 35;
    p.mu = 4;
    p.coderate = 2/3;
elseif MCS == 4
    snr = [20 25 30 35 40 42.5 45] %[-5:5:35]; % [0 10 20 30 50 80 100 200]; %0:5:35; % 0:10:50;
   % high_snr = 40;
    p.mu= 6;     % 2^mu-QAM.
    p.coderate = 3/4; % 1/3; Puncturing pattern is defined for 1/3 & 3/4 only
end

N = p.M*p.K;
noblk = 2;
Ncp = p.K;
N_tot = (N+Ncp);
N_2blk = (N+2*Ncp);
p.Delta_k = 2; %6;
Delta_k = p.Delta_k;
Fs = 1.92e6;
%path_delays = [0 1 3 5 6 11 37 53 115]+1; % [1 3 4 6 8 10 11 16 17];
%path_delays = [0 1 2 3 4 6 13 18 38]+1;
%AveGaindb =   [-1 -1 -1 0 0 0 -3 -5 -7];
path_delays = [0:8]+1;   %[0:5]+1;
numPaths = length(path_delays); %[10 15 20 25]; %[10 15 20 25 50]; 
AveGaindb = (linspace(0,-20, numPaths));
numPaths_delays = length(path_delays); %[10 15 20 25]; %[10 15 20 25 50]; 
numEndTaps = 0;
numTaps = numPaths_delays; %18+numEndTaps;  
fft_path_delays = [1:numTaps-numEndTaps (N-numEndTaps+1):N]; %[1:numTaps-10 N-9:N];% N-9:N];
fft_path_delays_cpblk = [1:numTaps-numEndTaps (N_2blk-numEndTaps+1):N_2blk]; %[1:numTaps-10 N_2blk-9:N_2blk];
%AveGaindb = (linspace(0,-20, numTaps));


p_d = p; 
p_d.Mp = 0; p_d.Md = 7;


p_ofdm = p;
p_ofdm.M = 1;
p_ofdm.K = p.M*p.K;
p_ofdm.Kon = p.M*p.Kon;
%p_ofdm.mu = p.mu;
p_ofdm.mode = 'ofdm';
p_ofdm.Delta_k = p.Delta_k*p.M;

p_ofdm.a = 0;
p_ofdm.L = p_ofdm.K;
p_ofdm.pulse = 'rc_td';

%chan_typ = 'time variant';
chan_typ = 'block fading';


nT = 2;     % no. of Tx antennas
nR = 2;     % no. of Rx antennas

% Doppler shift
%fd = 70; % 5Hz EPA, 70Hz EVA, 300Hz ETU
if fd == 0
    chan_typ = 'block fading';
else
    chan_typ = 'time variant';
end


%its = 100; % no. of simulation iterations

no_PIC_iterations = 3; 
no_PIC_iterations_Z = 5;
no_PIC_iterations_Zcp = 5;



%% Initial matrices and configurations

F = get_DFT_matrix(p);
F_M = dftmtx(p.M)/sqrt(p.M);
A = get_transmitter_matrix(p);
A_ofdm = get_transmitter_matrix(p_ofdm);

p_cp = p; p_cp.M = 1; p_cp.K = N+Ncp;
F_cp = get_DFT_matrix(p_cp);
F_N = get_DFT_matrix(p);
F_cpblk = dftmtx(Ncp+Ncp+N)/sqrt(Ncp+Ncp+N);

Fp = get_DFT_at_pilots('wider',p);
Fp_old = F_N; %get_DFT_at_pilots('exact',p);
Np = N; %N/Delta_k;
Fp_oldA = Fp_old*A;


if (Delta_k == 1)
   Fp = F; 
end
size_pfreqs = size(Fp,1);

    
% Power Delay Profile
P = zeros(numPaths,nT,nR);
for iR = 1:nR
    for iT = 1:nT
        
        a = 10.^(AveGaindb/10);
        %a = ones(1,numPaths);
        P(:,iT,iR) = a./sum(a); %a./sqrt(sum(a.^2));
    end
end
    
F_Lcpblk = F_cpblk(:,fft_path_delays_cpblk);

F_Lcp = F_cp(:,fft_path_delays);
Fbig_cp = kron(eye(nT),F_Lcp);
Fbig_Lcpblk = kron(eye(nT),F_Lcpblk);
Fbig_N = kron(eye(nT),F_N(:,fft_path_delays));
Fbig_Lcpblk_nTnR = kron(eye(nT*nR),F_Lcpblk);
F_L_N = F_N(:,fft_path_delays);


F_L = Fp(:,[fft_path_delays]);
Fbig = kron(eye(nT),F_L);

Z = kron(F_M,eye(p.K));
Z_mimo = kron(eye(nR),Z);

%F_L_old = Fp_old(:,[1:numPaths N-numPaths/2+1:end]);
F_L_old = Fp_old(:,fft_path_delays);
Fbig_old = kron(eye(nT),F_L_old);

I_cp = eye(N);
I_cp = [I_cp(end-Ncp+1:end,:); I_cp];
IcpA = blkdiag(A(end-Ncp+1:end,:), I_cp*A);
IcpA_ofdm = blkdiag(A_ofdm(end-Ncp+1:end,:), I_cp*A_ofdm);


conf.J_cp2 = [zeros(N,Ncp) eye(N) zeros(N,Ncp)];
conf.J = [zeros(N+Ncp,Ncp) eye(N+Ncp)];
conf.P = P;
conf.numPaths = numPaths;
conf.nT = nT; conf.nR = nR;
conf.path_delays = fft_path_delays;
conf.path_delays_cpblk = fft_path_delays_cpblk;

conf.Delta_k = Delta_k;
conf.noblk = noblk;
conf.N = N;
conf.Ncp = Ncp;
conf.N_2blk = N_2blk;
conf.N_tot = N_tot;
conf.F_L = F_Lcpblk;
conf.A = A;
conf.F_cpblk = F_cpblk;
conf.FcpIcpA = F_cpblk*IcpA;
conf.FA = F*A;
conf.Fbig_cp = Fbig_cp;
conf.F_Lcp = F_Lcp;
conf.Fbig_old = Fbig_old;
conf.F_L_old = F_L_old;
conf.F_L_N = F_L_N;
conf.JF = conf.J*F_cpblk'; 
conf.JF2 = conf.J_cp2*F_cpblk';
conf.FIcpA = F_cpblk*blkdiag(A(end-Ncp+1:end,:),I_cp*A);
conf.FIA = conf.FIcpA(:,end-N+1:end);
conf.F = F;
conf.M = p.M;
conf.K = p.K;
conf.Z = Z;

conf_ofdm = conf;
conf_ofdm.A = A_ofdm;
conf_ofdm.FcpIcpA = F_cpblk*IcpA_ofdm;
conf_ofdm.FA = F*A_ofdm;

JF = conf.JF; 
ZJF2 = Z*conf.JF2;

InT_cp = [];
for iT= 1:nT
    InT_cp = [InT_cp eye(N+2*Ncp)];
end
% preH_long = kron(eye(nR),InT_cp);

%% Pilot initializations

off_carriers = [1:k_start p.K-k_end+1:p.K];
Dd1_indices = ones(p.K,p.M);
Dd1_indices(off_carriers,:) = 0; 

Dd1_indices_ofdm = reshape(Dd1_indices.',N,1);
Dd1_indices = Dd1_indices(:);


Dd = zeros(p.K,p.M,nT); 
Dd_indx = zeros(p.K,p.M,nT); Dd1_indx = zeros(p.K,p.M,nT); 

Dd_indx_ofdm = zeros(p.K*p.M,1,nT); Dd1_indx_ofdm = zeros(p.K*p.M,1,nT);
Dp = zeros(p.K,p.M,nT);
%Dp1 = Dp; Dp2 = Dp;
xp_iT = zeros(N,nT); Xp_iT = zeros(N,nT);
xp_cp = zeros(N+Ncp,nT); Xp_cpiT = xp_cp;
xp_2blk = []; xp_ofdm_2blk = [];   
Xp = zeros(size_pfreqs,size_pfreqs*nT);
Xp_cp = zeros(N_tot,N_tot*nT);
Xp_old = zeros(Np,Np*nT);
xp_ofdm = zeros(N,nT); xp_cp_ofdm = zeros(N+Ncp,nT); Xp_ofdm = zeros(N,nT);
for blk = 1:noblk
    if blk == 1
        xp_cp = zeros(N+Ncp,nT);
        xp_cp_ofdm = zeros(N+Ncp,nT);
    else
    
        for iT = 1:nT
            [~, Dp(:,:,iT)] = do_pilot_symbols_time(p, Dd(:,:,iT),Delta_k,'Rand');
            %Dp([1:k_start end-k_end+1:end],:,iT) = 0;
            
            Dd_indx(:,:,iT) = abs((abs(Dp(:,:,iT)) > 0) - 1.0);
            Dd_indx([1:k_start end-k_end+1:end],:,iT) = 0;
            
            Dd1_indx(:,:,iT) = ones(p.K,p.M);
            Dd1_indx([1:k_start end-k_end+1:end],:,iT) = 0;                       
            
            Dd_indx_ofdm(:,:,iT) = reshape(Dd_indx(:,:,iT).', [N 1]);
            Dd1_indx_ofdm(:,:,iT) = reshape( Dd1_indx(:,:,iT).', [N 1]);            
            Dp_ofdm = reshape(Dp(:,:,iT).',N,1);     
             
            xp_iT(:,iT) =  do_modulate(p, Dp(:,:,iT));
            Xp_iT(:,iT) = fft_u(xp_iT(:,iT));
                
            xp_cp(:,iT) = [xp_iT(end-Ncp+1:end,iT); xp_iT(:,iT)];
            Xp_cpiT(:,iT) = F_cp*(xp_cp(:,iT));
            
            xp_ofdm(:,iT) = do_modulate(p_ofdm, Dp_ofdm);
            xp_cp_ofdm(:,iT) = [xp_ofdm(end-Ncp+1:end,iT); xp_ofdm(:,iT)];
            Xp_ofdm(:,iT) = fft_u(xp_ofdm(:,iT));
            
            % Xp in wide matrix form
            Xp(1:N,(N*(iT-1)+1):N*iT) = diag(Xp_iT(:,iT));
            Xp_cp(1:N_tot,(N_tot*(iT-1)+1):N_tot*iT) = diag(Xp_cpiT(:,iT));
            Xp_old(1:Np,(Np*(iT-1)+1):Np*iT) = diag(Fp_old*xp_iT(:,iT));
                               
        end
    end    
    xp_2blk = [xp_2blk; xp_cp];
    xp_ofdm_2blk = [xp_ofdm_2blk; xp_cp_ofdm];
end
% Dd_indx_on = Dd_indx(:,:,1);
% Dd_indx_on([1:k_start end-k_end+1:end],:) = [];
% 
% Dd1_indx_on = Dd1_indx(:,:,1);
% Dd1_indx_on([1:k_start end-k_end+1:end],:) = [];
% 
% Dd_indx_ofdm_on_iT = reshape( Dd_indx_on.', [numel(Dd_indx_on) 1]);
% Dd1_indx_ofdm_on_iT = reshape( Dd1_indx_on.', [numel(Dd1_indx_on) 1]);
% Dd_indx_on_ofdm = logical(repmat(Dd_indx_ofdm_on_iT, [nT 1]));
% Dd1_indx_on_ofdm = logical(repmat(Dd1_indx_ofdm_on_iT, [nT 1]));

x_test = xp_2blk(N_tot+Ncp+1:N_tot+N_tot,1);
pwr_x_test = abs(diag(x_test*x_test'));
numerator = 0; denomintr = 0;
for i=1:length(x_test)
    numerator = numerator + i*pwr_x_test(i);
    denomintr = denomintr + pwr_x_test(i);
end
wmean_pwr = numerator/denomintr;
n_H_gfdm = Ncp + floor(wmean_pwr);
n_H_ofdm = Ncp + N/2;
n_H_gfdm_cp = N_tot/2;

n_H = N_tot + [n_H_gfdm n_H_ofdm n_H_gfdm_cp];


XpF = Xp*Fbig_N;   
XpF_old = Xp_old*Fbig_old;  
if nT == 1
    xp_NcpNcpN = xp_2blk(end-N-2*Ncp+1:end);
    Xp_2blk = diag(F_cpblk*xp_NcpNcpN);
    
    xp_ofdm_NcpNcpN = xp_ofdm_2blk(end-N_tot-Ncp+1:end);
    Xp_2blk_ofdm = diag(F_cpblk*xp_ofdm_NcpNcpN);
        
else    
    % Reshape the pilots signal into wide matrix of diagonals
    wXp_2blk = fft_u(xp_2blk(end-N_tot-Ncp+1:end,:));
    dXp_2blk = diag(wXp_2blk(:));
    Xp_2blk = InT_cp*dXp_2blk;
    
    wXp_2blk_ofdm = fft_u(xp_ofdm_2blk(end-N_tot-Ncp+1:end,:));
    dXp_2blk_ofdm = diag(wXp_2blk_ofdm(:));
    Xp_2blk_ofdm = InT_cp*dXp_2blk_ofdm;
end
XpF_cp = JF*Xp_2blk*Fbig_Lcpblk;
XpF_cp_ofdm = JF*Xp_2blk_ofdm*Fbig_Lcpblk;
R_PsiPsi_Z = 0; %calc_Rpsipsi('R_PsiPsi_Z', nT, nR, conf, Dd_indx, 0);
R_PsiPsi_cp = calc_Rpsipsi('R_PsiPsi_CP', nT, nR, conf, Dd_indx, Dd1_indices);
R_PsiPsi_ofdm = calc_Rpsipsi('R_PsiPsi_CP', nT, nR, conf_ofdm, Dd_indx_ofdm, Dd1_indices_ofdm);

XpF_cp_2_Zm = zeros(p.K,nT*numPaths,p.M);
XpF_old_Zm = zeros(N,nT*numPaths,p.M); 
XpF_old_Zm_white = zeros(N,nT*numPaths,p.M); 


for m = 1:p.M
    Z_m = kron(F_M(m,:),eye(p.K));
    XpF_cp_2_Zm(:,:,m) = Z_m*conf.JF2*Xp_2blk*Fbig_Lcpblk;
    XpF_old_Zm(:,:,m) = F'*XpF_old;  
end

Dd_i = Dd_indx(:,:,1);
Dd1_i = Dd1_indx(:,:,1);
sigma2d = (abs( Dd_i(:) ) > 0) + 0;
sigma2d_blk1 = (abs( Dd1_i(:) ) > 0) + 0;
R_dd = kron(eye(nT),diag(sigma2d));
R_dd_blk1 = kron(eye(nT),diag(sigma2d_blk1));
FA_Rdd_AF_tx1 = F*A*diag(sigma2d)*A'*F';

Dd_i_ofdm = Dd_indx_ofdm(:,:,1);
Dd1_i_ofdm = Dd1_indx_ofdm(:,:,1);
sigma2d_ofdm = (abs( Dd_i_ofdm(:) ) > 0) + 0;
sigma2d_blk1_ofdm = (abs( Dd1_i_ofdm(:) ) > 0) + 0;
R_dd_ofdm = kron(eye(nT),diag(sigma2d_ofdm));
R_dd_blk1_ofdm = kron(eye(nT),diag(sigma2d_blk1_ofdm));

[on_bins, on_bins1, off_bins, on_bins_single] = get_on_bins(Xp_iT, N, nT,nR);
[on_bins_ofdm, on_bins1_ofdm, off_bins_ofdm, on_bins_single_ofdm] = get_on_bins(Xp_ofdm, N, nT,nR);


pilots_conf.off_bins = off_bins;
pilots_conf.on_bins = on_bins;
pilots_conf.on_bins1 = on_bins1;
pilots_conf.on_bins_single = on_bins_single;
pilots_conf.FA_Rdd_AF_tx1 = FA_Rdd_AF_tx1;
pilots_conf.XpF = XpF;
pilots_conf.Xp = Xp;
pilots_conf.R_dd = R_dd;
pilots_conf.R_dd_blk1 = R_dd_blk1;
pilots_conf.XpF_cp = XpF_cp;
pilots_conf.XpF_old_Zm= XpF_old_Zm;
pilots_conf.XpF_cp_2_Zm = XpF_cp_2_Zm;
pilots_conf.Xp_iT = Xp_iT;
pilots_conf.Dp = Dp;
pilots_conf.Dd_indx = Dd_indx;
pilots_conf.Dd1_indx = Dd1_indx;

pilots_conf_ofdm.off_bins = off_bins_ofdm;
pilots_conf_ofdm.on_bins = on_bins_ofdm;
pilots_conf_ofdm.on_bins1 = on_bins1_ofdm;
pilots_conf_ofdm.on_bins_single = on_bins_single_ofdm;
pilots_conf_ofdm.Xp_iT = Xp_ofdm;
pilots_conf_ofdm.R_dd = R_dd_ofdm;
pilots_conf_ofdm.R_dd_blk1 = R_dd_blk1_ofdm;
pilots_conf_ofdm.FA_Rdd_AF_tx1 = zeros(size(FA_Rdd_AF_tx1));
pilots_conf_ofdm.XpF = zeros(size(XpF));
pilots_conf_ofdm.Dd_indx = Dd_indx_ofdm;
pilots_conf_ofdm.Dd1_indx = Dd1_indx_ofdm;
% pilots_conf_ofdm.Dd1_indx_on = Dd1_indx_on_ofdm;
% pilots_conf_ofdm.Dd_indx_on = Dd_indx_on_ofdm;


p_blk1 = p;
p_blk1.Mp = 0; p_blk1.Md = p.M;   
%p.Mp = 2; p.Md = 5;

%% Simulation

% Initialize turbo coding
Nd_blk1 = length(find(sigma2d_blk1));  % Number of QAM symbols per antenna
[PCCC_blk1, interleaver_blk1, bicm_blk1] = find_interleaver_size(p, Nd_blk1, nT); % PCCC_blk1.code_param.tail_pattern
Nd_blk2 = length(find(sigma2d)); % Number of QAM symbols per antenna
[PCCC_blk2, interleaver_blk2, bicm_blk2] = find_interleaver_size(p, Nd_blk2, nT);

no_estimattion_coeff = numPaths*nT*nR;
no_observations = length(find(abs(Dp(:,1,:))>1e-8));
if no_observations < no_estimattion_coeff
   warning('Number of estimation coefficients are larger than the number of observations \n no_est_coeff: %d, no_obs: %d', no_estimattion_coeff,no_observations) 
else
   fprintf('no_est_coeff: %d, no_obs: %d \n',no_estimattion_coeff,no_observations)
end
if ~isfield(conf,'test_snr')
   warning('off','all')
end
rng(rnd_seed);
for si=1:length(snr)
        
    Cstamp = clock;
    texttodisp = sprintf('L = %d, SNR = %.3f, on %02d.%02d @ %02d:%02d started with %d iterations',numPaths,snr(si),Cstamp(2),Cstamp(3),Cstamp(4),Cstamp(5),its);
    disp(texttodisp)       
    
        
   % p.Mp = 2; p.Md = 5;
    for i=1:its
        
        % Initializations        
        %s = zeros(N-p.Mp*p.K/Delta_k,nT); dd = zeros(N-p.Mp*p.K/Delta_k,nT);
        Dd = zeros(p.K,p.M,nT); %Dp = zeros(p.K,p.M,nT); 
        D = zeros(p.K,p.M,nT);  Dd_ofdm = zeros(p.K*p.M,1,nT); Dd1_ofdm = zeros(p.K*p.M,1,nT);
        
        xd_iT = zeros(N,nT); Xd_iT = zeros(size_pfreqs,nT);
        
        Xd = zeros(size_pfreqs,size_pfreqs*nT);
        X = zeros(size_pfreqs,size_pfreqs*nT);
        x_iT = zeros(N,nT); X_iT = zeros(size_pfreqs,nT);   
        x_iT_ofdm = zeros(N,nT); X_iT_ofdm = zeros(size_pfreqs,nT);        
        
        x_cp = zeros(N+Ncp,nT); X_cpiT = x_cp;
        
        x_cp_ofdm = zeros(N+Ncp,nT);
        
        x_ofdm_2blk = []; 
        x_2blk = [];
        xp_2blk = []; xp_ofdm_2blk = [];     
        for blk = 1:noblk
            % Signal generation for Tx antennas
            if blk == 1               
                
                dd = zeros(Nd_blk1,nT);
                Dp0 = zeros(size(Dp));                
                
                % Random interleavers
                code_interleaver_blk1 = randperm(interleaver_blk1)-1;
                bicm_interleaver_blk1 = randperm(bicm_blk1)-1;
                
                % Random binary bits
                b_blk1 = round(rand([1 interleaver_blk1]));
                
                % Encoder
                [dd_nT, bc_nT_blk1 ] = do_encode_pccc(b_blk1, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
                
                dd_nT1 = dd_nT;
                         
            elseif blk == 2
                
                
                dd = zeros(Nd_blk2,nT);
                Dp0 = Dp;         
                
                % Random interleavers
                code_interleaver_blk2 = randperm(interleaver_blk2)-1;
                bicm_interleaver_blk2 = randperm(bicm_blk2)-1;
                
                % Random binary bits
                b_blk2 = round(rand([1 interleaver_blk2]));
                
                % Encoder
                [dd_nT, bc_nT_blk2 ] = do_encode_pccc(b_blk2, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
                
                dd_nT2 = dd_nT;
            end
            
            
            dd = reshape(dd_nT,size(dd));
            
            
            for iT = 1:nT
                if blk == 1
                   
                    Dd(:,:,iT) = do_map_p(p_blk1, dd(:,iT)); % map them to the D matrix
                         
                    Dp0(off_carriers,:,iT) = 0;
                    
                    Dd1 = Dd;
                    Dd1_ofdm(:,:,iT) = reshape(Dd1(:,:,iT).',N,1);
                else                                       
                
                    Dd(:,:,iT) = do_map_p(p, dd(:,iT)); % map them to the D matrix                    
                            
                    Dp0(off_carriers,:,iT) = 0;
                end        
                
                % Pilots and Data signal
                D(:,:,iT) =  Dd(:,:,iT) + Dp0(:,:,iT);
                x_iT(:,iT) = do_modulate(p, D(:,:,iT));
                %X_iT(:,iT) = Fp*(x_iT(:,iT));
                
                x_cp(:,iT) = [x_iT(end-Ncp+1:end,iT); x_iT(:,iT)];
                %X_cpiT(:,iT) = F_cp*(x_cp(:,iT));
                
%                 xd_iT(:,iT) =  do_modulate(p, Dd(:,:,iT));
%                 Xd_iT(:,iT) = fft_u(xd_iT(:,iT));
                
                % OFDM
                D_ofdm = reshape(D(:,:,iT).',N,1);                
                % Dd_ofdm(:,:,iT) = reshape(Dd(:,:,iT).',N,1);
                
                x_iT_ofdm(:,iT) = do_modulate(p_ofdm, D_ofdm);
                x_cp_ofdm(:,iT) = [x_iT_ofdm(end-Ncp+1:end,iT); x_iT_ofdm(:,iT)];
                                
                % X(1:size_pfreqs,(size_pfreqs*(iT-1)+1):size_pfreqs*iT) = diag(X_iT(:,iT));
                % Xd(1:size_pfreqs,(size_pfreqs*(iT-1)+1):size_pfreqs*iT) = diag(Xd_iT(:,iT));
                
                
            end
            x_2blk = [x_2blk; x_cp];            
            x_ofdm_2blk = [x_ofdm_2blk; x_cp_ofdm];
        end             
        
                        
        % Random channel  
        h = get_mimo_channel(P,Fs,numPaths,path_delays,nT,nR,chan_typ,fd);                   


        % Apply the MIMO Channel to the multiple GFDM data blocks
        % Note: CP is necessary to create a circulant channel
        if strcmp(chan_typ,'time variant')
            xch_2blk = do_mimo_channel('Filter time variant',x_2blk, h, numPaths, nT, nR, N_tot, noblk);

            % Get the Genie channel (Must be after channel filtering)
            [h_gains, h_gains_m] = get_genie_channel(h,chan_typ,path_delays,numTaps,nT,nR,N,n_H);
            vH_genie_gfdm = reshape(sqrt(N)*Fbig_N*h_gains{1},[(N)*nT*nR 1]);
            %vH_genie_ofdm = reshape(sqrt(N)*Fbig_N*h_gains{2},[(N)*nT*nR 1]);
            vH_genie_gfdm_cp = reshape(sqrt(N)*Fbig_N*h_gains{3},[(N)*nT*nR 1]);
            vH_genie_m = reshape(sqrt(N)*Fbig_N*h_gains_m,[(N)*nT*nR 1]);

            xch_ofdm_2blk = do_mimo_channel('Filter time variant',x_ofdm_2blk, h, numPaths, nT, nR, N_tot, noblk);

            [h_gains, h_gains_m] = get_genie_channel(h,chan_typ,path_delays,numTaps,nT,nR,N,n_H);
            %vH_genie_gfdm = reshape(sqrt(N)*Fbig_N*h_gains{1},[(N)*nT*nR 1]);
            vH_genie_ofdm = reshape(sqrt(N)*Fbig_N*h_gains{2},[(N)*nT*nR 1]);
            %vH_genie_gfdm_cp = reshape(sqrt(N)*Fbig_N*h_gains{3},[(N)*nT*nR 1]);
            vH_genie_m_ofdm = reshape(sqrt(N)*Fbig_N*h_gains_m,[(N)*nT*nR 1]);
        else
            xch_2blk = do_mimo_channel('Filter',x_2blk, h_gains, numPaths, nT, nR, N_tot, noblk);         
            xch_ofdm_2blk = do_mimo_channel('Filter',x_ofdm_2blk, h_gains, numPaths, nT, nR, N_tot, noblk);        
        end
        
        % Get the Genie channel (Must be after channel filtering)
        % [h_gains, h_gains_m] = get_genie_channel(h,chan_typ,path_delays,numTaps,nT,nR,N);
        %vH_genie = reshape(sqrt(N)*Fbig_N*h_gains,[(N)*nT*nR 1]);
        %vH_genie_m = reshape(sqrt(N)*Fbig_N*h_gains_m,[(N)*nT*nR 1]);
        
        if isfield(conf,'test_snr')
            snr(si) = conf.test_snr;
            warning('Test SNR is %d',conf.test_snr);
        end
    
        % AWGN
        w_2blk = awgn(zeros(noblk*N_tot,nR),snr(si));
        
        y_2blk = xch_2blk + w_2blk; % Two blocks GFDM
        y_ofdm_2blk = xch_ofdm_2blk + w_2blk; % Two blocks OFDM
                        
        y = y_2blk(end-N+1:end,:); % 2nd block GFDM
        y1 = y_2blk(Ncp+1:Ncp+N,:); % 1st block GFDM
       
        y_ofdm = y_ofdm_2blk(end-N+1:end,:); % 2nd block OFDM
        y1_ofdm = y_ofdm_2blk(Ncp+1:Ncp+N,:); % 1st block OFDM
        
        ycp = y_2blk(end-N_tot+1:end,:); % 2nd block with CP, GFDM
        ycp_ofdm = y_ofdm_2blk(end-N_tot+1:end,:); % 2nd block with CP, OFDM
                
        % Receive signal at pilot subcarriers
        Y = Fp*(y);   
        
        Y2_ofdm = fft_u(y_ofdm);
        Y1_ofdm = fft_u(y1_ofdm);
        
        Y_all = fft_u(y);
        Y1_all = fft_u(y1);
        
        Y_all_v = [Y1_all(:) Y_all(:)];         

        [R_PsiPsi, ~, ~] = calc_ce_mse_mimo('R_PsiPsi', Y, snr(si), N, nT, nR, Fbig_old, F_L_old, Fp_oldA, XpF_old, P, numPaths, Dd, Delta_k, A); 
        R_PsiPsi(abs(R_PsiPsi)<1e-6) = 0;              

        
 %%      GFDM LMMSE ESTIMATION        
         [h_lmmse, R_hh_err_lmmse] = do_Chan_Estimation_mimo('LMMSE', Fp_old*y, snr(si), N, nT, nR, Fbig_old, F_L_old, Fp_oldA, XpF_old, P, numPaths, path_delays, Dd, Delta_k, ((R_PsiPsi)), 0, p.M, 0, on_bins1);         
         wh_lmmse= reshape(h_lmmse,[numPaths nT*nR]);
         wH_lmmse = sqrt(N)*F_L_N*(wh_lmmse);
         vH_lmmse = wH_lmmse(:);          
         
 %%      OFDM CP LMMSE        
         [vH_lmmse_cp_ofdm, R_hh_err_ofdm] = do_Chan_Estimation_mimo('LMMSE_CP', ycp_ofdm, snr(si), N_2blk, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp_ofdm, P, numPaths, fft_path_delays, Dd, Delta_k, R_PsiPsi_ofdm,Ncp);
                 
        
         
         
 %%      GFDM CP LMMSE        
         [vH_lmmse_cp, R_hh_err_cp] = do_Chan_Estimation_mimo('LMMSE_CP', ycp, snr(si), N_2blk, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp, P, numPaths, fft_path_delays, Dd, Delta_k, R_PsiPsi_cp,Ncp);
         vH_lmmse_in = vH_lmmse_cp;   
         
         [strucH_lmmse_Zm, R_hh_lmmsez] = do_Chan_Estimation_mimo('LMMSE_Zm', y, min(snr(si), max_ce_input_snr), N, nT, nR, Fbig_old, F_L_old, 0, XpF_old_Zm, P, numPaths, fft_path_delays, Dp, Delta_k, R_PsiPsi_Z,Ncp,p.M,p.K);
         vH_lmmse_Zm = strucH_lmmse_Zm.vH_lmmse;
         
         [vH_lmmse_cp2_m, R_hh_lmmsezcp] = do_Chan_Estimation_mimo('LMMSE_CPZm', ycp(1:N,:), min(snr(si), max_ce_input_snr), N+Ncp, nT, nR, Fbig_cp, F_Lcp, 0, XpF_cp_2_Zm, P, numPaths, fft_path_delays, Dp, Delta_k, strucH_lmmse_Zm,Ncp,p.M,p.K);
         vH_lmmse_cpZm = vH_lmmse_cp2_m; 
         
%          mean(abs(vH_lmmse_Zm(on_bins) - vH_big(on_bins) ).^2)

%        % Sequential LMMSE with whitenning
%
%          if white_LMMSE_Zm
%              VarN = 1/(10^(snr(si)/10));
%              whitener = sqrtm( R_PsiPsi(1:N,1:N) + VarN*eye(N) );
%              for m = 1:p.M
%                  XpF_old_Zm_white(:,:,m) = sqrt(N)*ifft((sparse(whitener'))\XpF_old);
%              end
%              y_white = sqrt(N)*ifft(((whitener'))\fft_u(y));
%              snr_white = 0;
%              
%              [strucH_lmmse_Zm, R_hh_lmmsez] = do_Chan_Estimation_mimo('LMMSE_Zm', y_white, snr_white, N, nT, nR, Fbig_old, F_L_old, 0, XpF_old_Zm_white, P, numPaths, fft_path_delays, Dp, Delta_k, R_PsiPsi_Z,Ncp,p.M,p.K);
%              vH_lmmse_Zm_white = strucH_lmmse_Zm.vH_lmmse;
%              
% %              mean(abs(vH_lmmse_Zm_white(on_bins) - vH_big(on_bins) ).^2)
%          end     
         d_genie = [reshape(Dd1,[N 1 nT]); reshape(Dd,[N 1 nT])];
         [vH_lmmse_3, dhat_gfdm_lmmse_cp, Heq_3 ] = do_CE_MMSE_PIC('LMMSE_CP',p, conf, pilots_conf, no_PIC_iterations, snr(si), vH_lmmse_in, R_hh_err_cp, Dd1_indices, ycp, Y_all_v, d_genie(:), vH_genie_gfdm_cp );
         
         [vH_lmmse_3z, dhat_gfdmz, Heq_3z ] = do_CE_MMSE_PIC('LMMSE_Zm',p, conf, pilots_conf, no_PIC_iterations_Z, [max_ce_input_snr snr(si)], vH_lmmse_Zm, R_hh_lmmsez, Dd1_indices, ycp(Ncp+1:end,:), Y_all_v, Dd(:), vH_genie_gfdm );              
         
         
         [vH_lmmse_3z_cp, dhat_gfdmzcp, Heq_3z_cp ] = do_CE_MMSE_PIC('LMMSE_CPZm',p, conf, pilots_conf, no_PIC_iterations_Zcp, [max_ce_input_snr snr(si)], vH_lmmse_cpZm, R_hh_lmmsezcp, Dd1_indices, ycp, Y_all_v, d_genie(:), vH_genie_gfdm_cp );

         
%% Symbol Detection
if detection_flag
    
         % OFDM
         i_blk = 1; ofdm_flag = true; % R_hh_err_ofdm = zeros(numPaths*nT*nR,numPaths*nT*nR);  
         [d_hat_soft_ofdm, Heq_ofdm] = do_detect_soft_sym(i_blk, p, conf_ofdm, pilots_conf_ofdm, Y1_ofdm, snr(si), vH_lmmse_cp_ofdm, R_hh_err_ofdm, N, nT, nR, ofdm_flag);
         i_blk = 2;
         [d_hat_soft_ofdm, Heq_ofdm] = do_detect_soft_sym(i_blk, p, conf_ofdm, pilots_conf_ofdm, Y2_ofdm, snr(si), vH_lmmse_cp_ofdm, R_hh_err_ofdm, N, nT, nR, ofdm_flag, d_hat_soft_ofdm, Heq_ofdm);

         % GFDM, LMMSE Channel Estimation 
         i_blk = 1; ofdm_flag = false; 
         [d_hat_soft_lmmse, Heq_lmmse] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, snr(si), vH_lmmse, R_hh_err_lmmse, N, nT, nR, ofdm_flag);
         i_blk = 2; 
         [d_hat_soft_lmmse, Heq_lmmse] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, snr(si), vH_lmmse, R_hh_err_lmmse, N, nT, nR, ofdm_flag, d_hat_soft_lmmse, Heq_lmmse);
         
         % GFDM, CP-LMMSE Channel Estimation 
         i_blk = 1; ofdm_flag = false; 
         [d_hat_soft_cplmmse, Heq_cplmmse] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, snr(si), vH_lmmse_cp, R_hh_err_cp, N, nT, nR, ofdm_flag);
         i_blk = 2; 
         [d_hat_soft_cplmmse, Heq_cplmmse] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, snr(si), vH_lmmse_cp, R_hh_err_cp, N, nT, nR, ofdm_flag, d_hat_soft_cplmmse, Heq_cplmmse);
         
         % GFDM, LMMSE-Z Channel Estimation 
         i_blk = 1; ofdm_flag = false; 
         [d_hat_soft_lmmsez, Heq_lmmsez] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, min(snr(si), max_ce_input_snr), vH_lmmse_Zm, R_hh_lmmsez, N, nT, nR, ofdm_flag);
         i_blk = 2; 
         [d_hat_soft_lmmsez, Heq_lmmsez] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, min(snr(si), max_ce_input_snr), vH_lmmse_Zm, R_hh_lmmsez, N, nT, nR, ofdm_flag, d_hat_soft_lmmsez, Heq_lmmsez);
         
         
         % GFDM, CP-LMMSE-Z Channel Estimation 
         i_blk = 1; ofdm_flag = false; 
         [d_hat_soft_cplmmsez, Heq_cplmmsez] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, min(snr(si), max_ce_input_snr), vH_lmmse_cpZm, R_hh_lmmsezcp, N, nT, nR, ofdm_flag);
         i_blk = 2; 
         [d_hat_soft_cplmmsez, Heq_cplmmsez] = do_detect_soft_sym(i_blk, p, conf, pilots_conf, Y_all_v, min(snr(si), max_ce_input_snr), vH_lmmse_cpZm, R_hh_lmmsezcp, N, nT, nR, ofdm_flag, d_hat_soft_cplmmsez, Heq_cplmmsez);


%% Decode
         
         EsNo = 1;%10^(snr(si)/10);
         
         % GFDM LMMSE
         [bh_lmmse1, bit_llr_lmmse1] = do_decode_pccc( d_hat_soft_lmmse{1}, b_blk1, EsNo, Heq_lmmse{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_lmmse2, bit_llr_lmmse2] = do_decode_pccc( d_hat_soft_lmmse{2}, b_blk2, EsNo, Heq_lmmse{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % GFDM LMMSE-Z
         [bh_lmmsez1, bit_llr_lmmsez1] = do_decode_pccc( d_hat_soft_lmmsez{1}, b_blk1, EsNo, Heq_lmmsez{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_lmmsez2, bit_llr_lmmsez2] = do_decode_pccc( d_hat_soft_lmmsez{2}, b_blk2, EsNo, Heq_lmmsez{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % GFDM CP-LMMSE-Z
         [bh_cplmmsez1, bit_llr_cplmmsez1] = do_decode_pccc( d_hat_soft_cplmmsez{1}, b_blk1, EsNo, Heq_cplmmsez{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_cplmmsez2, bit_llr_cplmmsez2] = do_decode_pccc( d_hat_soft_cplmmsez{2}, b_blk2, EsNo, Heq_cplmmsez{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % GFDM CP-LMMSE
         [bh_cplmmse1, bit_llr_cplmmse1] = do_decode_pccc( d_hat_soft_cplmmse{1}, b_blk1, EsNo, Heq_cplmmse{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_cplmmse2, bit_llr_cplmmse2] = do_decode_pccc( d_hat_soft_cplmmse{2}, b_blk2, EsNo, Heq_cplmmse{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);         
         
         % GFDM CP-LMMSE-PIC
         [bh_cepic1_cp, bit_llr_cepic1_cp] = do_decode_pccc( dhat_gfdm_lmmse_cp{1}, b_blk1, EsNo, Heq_3{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_cepic2_cp, bit_llr_cepic2_cp] = do_decode_pccc( dhat_gfdm_lmmse_cp{2}, b_blk2, EsNo, Heq_3{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % GFDM LMMSE-Z-PIC
         [bh_cepicz1, bit_llr_cepicz1] = do_decode_pccc( dhat_gfdmz{1}, b_blk1, EsNo, Heq_3z{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_cepicz2, bit_llr_cepicz2] = do_decode_pccc( dhat_gfdmz{2}, b_blk2, EsNo, Heq_3z{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % GFDM CP-LMMSE-Z-PIC
         [bh_cepicz1_cp, bit_llr_cepicz1_cp] = do_decode_pccc( dhat_gfdmzcp{1}, b_blk1, EsNo, Heq_3z_cp{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_cepicz2_cp, bit_llr_cepicz2_cp] = do_decode_pccc( dhat_gfdmzcp{2}, b_blk2, EsNo, Heq_3z_cp{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);
         
         % OFDM
         EsNo_ofdm = 1;
         [bh_ofdm1, bit_llr_ofdm1] = do_decode_pccc( d_hat_soft_ofdm{1}, b_blk1, EsNo_ofdm, Heq_ofdm{1}, PCCC_blk1, code_interleaver_blk1, bicm_interleaver_blk1);
         [bh_ofdm2, bit_llr_ofdm2] = do_decode_pccc( d_hat_soft_ofdm{2}, b_blk2, EsNo_ofdm, Heq_ofdm{2}, PCCC_blk2, code_interleaver_blk2, bicm_interleaver_blk2);

         
         % Coded-BER 
         bitErr_lmmse(si,i) = mean(( [bh_lmmse1 bh_lmmse2] ~= [b_blk1 b_blk2]));        
         bitErr_lmmsez(si,i) = mean(( [bh_lmmsez1 bh_lmmsez2] ~= [b_blk1 b_blk2]));        
         bitErr_cplmmsez(si,i) = mean(( [bh_cplmmsez1 bh_cplmmsez2] ~= [b_blk1 b_blk2]));        
         bitErr_cplmmse(si,i) = mean(( [bh_cplmmse1 bh_cplmmse2] ~= [b_blk1 b_blk2]));        
         bitErr_cepic_cp(si,i) = mean(( [bh_cepic1_cp bh_cepic2_cp] ~= [b_blk1 b_blk2]));        
         bitErr_cepicz(si,i) = mean(( [bh_cepicz1 bh_cepicz2] ~= [b_blk1 b_blk2]));         
         bitErr_cepicz_cp(si,i) = mean(( [bh_cepicz1_cp bh_cepicz2_cp] ~= [b_blk1 b_blk2])) ;   
         bitErr_ofdm(si,i) = mean(( [bh_ofdm1 bh_ofdm2] ~=  [b_blk1 b_blk2]   ));
         
	 % Coded-FER 
         fer_lmmse(si,i) = mean([~isequal(bh_lmmse1,b_blk1) ~isequal(bh_lmmse2,b_blk2)]);         
         fer_lmmsez(si,i) = mean([~isequal(bh_lmmsez1,b_blk1) ~isequal(bh_lmmsez2,b_blk2)]);              
         fer_cplmmsez(si,i) = mean([~isequal(bh_cplmmsez1,b_blk1) ~isequal(bh_cplmmsez2,b_blk2)]); 
         fer_cplmmse(si,i) = mean([~isequal(bh_cplmmse1,b_blk1) ~isequal(bh_cplmmse2,b_blk2)]);     
         fer_cepic_cp(si,i) = mean([~isequal(bh_cepic1_cp,b_blk1) ~isequal(bh_cepic2_cp,b_blk2)]);     
         fer_cepicz(si,i) = mean([~isequal(bh_cepicz1,b_blk1) ~isequal(bh_cepicz2,b_blk2)]);        
         fer_cepicz_cp(si,i) = mean([~isequal(bh_cepicz1_cp,b_blk1) ~isequal(bh_cepicz2_cp,b_blk2)]);   
         fer_ofdm(si,i) = mean([~isequal(bh_ofdm1,b_blk1) ~isequal(bh_ofdm2,b_blk2)]); 
         
         % Mutual information
         mi_lmmse(si,i) = mean([mutualinfo(bit_llr_lmmse1,bc_nT_blk1) mutualinfo(bit_llr_lmmse2,bc_nT_blk2)]);
         mi_lmmsez(si,i) = mean([mutualinfo(bit_llr_lmmsez1,bc_nT_blk1) mutualinfo(bit_llr_lmmsez2,bc_nT_blk2)]);          
         mi_cplmmsez(si,i) = mean([mutualinfo(bit_llr_cplmmsez1,bc_nT_blk1) mutualinfo(bit_llr_cplmmsez2,bc_nT_blk2)]); 
         mi_cplmmse(si,i) = mean([mutualinfo(bit_llr_cplmmse1,bc_nT_blk1) mutualinfo(bit_llr_cplmmse2,bc_nT_blk2)]);     
         mi_cepic_cp(si,i) = mean([mutualinfo(bit_llr_cepic1_cp,bc_nT_blk1) mutualinfo(bit_llr_cepic2_cp,bc_nT_blk2)]);   
         mi_cepicz(si,i) = mean([mutualinfo(bit_llr_cepicz1,bc_nT_blk1) mutualinfo(bit_llr_cepicz2,bc_nT_blk2)]);       
         mi_cepicz_cp(si,i) = mean([mutualinfo(bit_llr_cepicz1_cp,bc_nT_blk1) mutualinfo(bit_llr_cepicz2_cp,bc_nT_blk2)]);  
         mi_ofdm(si,i) = mean([mutualinfo(bit_llr_ofdm1,bc_nT_blk1) mutualinfo(bit_llr_ofdm2,bc_nT_blk2)]); 
        
         % Uncoded-SER
         symbols1 = do_qamdemodulate(dd_nT1,p.mu);
         symbols2 = do_qamdemodulate(dd_nT2,p.mu);        
         
         ser_ofdm(si,i) = calc_uncoded_ser(symbols1, symbols2, d_hat_soft_ofdm, Heq_ofdm, p);
         ser_lmmse(si,i) = calc_uncoded_ser(symbols1, symbols2, d_hat_soft_lmmse, Heq_lmmse, p);
         ser_lmmsez(si,i) = calc_uncoded_ser(symbols1, symbols2, d_hat_soft_lmmsez, Heq_lmmsez, p);
         ser_cplmmse(si,i) = calc_uncoded_ser(symbols1, symbols2, d_hat_soft_cplmmse, Heq_cplmmse, p);
         ser_cplmmsez(si,i) = calc_uncoded_ser(symbols1, symbols2, d_hat_soft_cplmmsez, Heq_cplmmsez, p);
         ser_cepicz(si,i) = calc_uncoded_ser(symbols1, symbols2, dhat_gfdmz, Heq_3z, p);
         ser_cepicz_cp(si,i) = calc_uncoded_ser(symbols1, symbols2, dhat_gfdmzcp, Heq_3z_cp, p);          
         ser_cepic_cp(si,i) = calc_uncoded_ser(symbols1, symbols2, dhat_gfdm_lmmse_cp, Heq_3, p);
         
         
	 	 
end      
   
%% Mean Square Error of Channel Estimations
	 
	 
         % Genie-Aided Channel and Data 
         mse_pic_new_genieHD(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3(on_bins,end) ).^2 );
         mse_pic_new_genieHD_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3(on_bins,end) ).^2 );         
         
         % Genie-Aided Data with imperfect channel
         mse_pic_new_genieD(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3(on_bins,end-1) ).^2 );
         mse_pic_new_genieD_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3(on_bins,end-1) ).^2 );    
         
         % Genie-Aided H&D, LMMSE-PIC-Z
         mse_pic_new_1z_genie(si,i) = mean( abs( vH_genie_gfdm(on_bins) - vH_lmmse_3z(on_bins,end) ).^2 );
         mse_pic_new_1z_genie_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3z(on_bins,end) ).^2 );
         
         % Genie-Aided H&D, CP-LMMSE-Z-PIC
         mse_pic_new_1zcp_genie(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3z_cp(on_bins,end) ).^2 );
         mse_pic_new_1zcp_genie_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3z_cp(on_bins,end) ).^2 );
         
         % CP-LMMSE-PIC
         mse_pic_new_1(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3(on_bins,no_PIC_iterations+1) ).^2 );
         mse_pic_new_1_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3(on_bins,no_PIC_iterations+1) ).^2 );
         
         % LMMSE-PIC-Z
         mse_pic_new_1z(si,i) = mean( abs( vH_genie_gfdm(on_bins) - vH_lmmse_3z(on_bins,no_PIC_iterations_Z+1) ).^2 );
         mse_pic_new_1z_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3z(on_bins,no_PIC_iterations_Z+1) ).^2 );
         
         % CP-LMMSE-Z-PIC
         mse_pic_new_1zcp(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3z_cp(on_bins,no_PIC_iterations_Zcp+1) ).^2 );
         mse_pic_new_1zcp_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3z_cp(on_bins,no_PIC_iterations_Zcp+1) ).^2 );
                  
         mse_pic_new_1zcp_genie(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_3z_cp(on_bins,end) ).^2 );
         mse_pic_new_1zcp_genie_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_3z_cp(on_bins,end) ).^2 );
                  
         mse_lmmse(si,i) = mean(abs(vH_lmmse(on_bins) - vH_genie_gfdm(on_bins) ).^2);                 
         mse_lmmse_m(si,i) = mean(abs(vH_lmmse(on_bins) - vH_genie_m(on_bins) ).^2);                 
         
         mse_lmmsecp(si,i) = mean( abs( vH_genie_gfdm_cp(on_bins) - vH_lmmse_cp(on_bins) ).^2 );
         mse_lmmsecp_m(si,i) = mean( abs( vH_genie_m(on_bins) - vH_lmmse_cp(on_bins) ).^2 );
         
         mse_lmmsecp_ofdm(si,i) = mean( abs( vH_genie_ofdm(on_bins) - vH_lmmse_cp_ofdm(on_bins) ).^2 );
         mse_lmmsecp_ofdm_m(si,i) = mean( abs( vH_genie_m_ofdm(on_bins) - vH_lmmse_cp_ofdm(on_bins) ).^2 );
        
         mse_lmmsez(si,i) = mean(abs(vH_lmmse_Zm(on_bins) - vH_genie_gfdm(on_bins) ).^2);                 
         mse_lmmsez_m(si,i) = mean(abs(vH_lmmse_Zm(on_bins) - vH_genie_m(on_bins) ).^2);                 
         
         mse_lmmsecpz(si,i) = mean(abs(vH_lmmse_cpZm(on_bins) - vH_genie_gfdm_cp(on_bins) ).^2);                 
         mse_lmmsecpz_m(si,i) = mean(abs(vH_lmmse_cpZm(on_bins) - vH_genie_m(on_bins) ).^2);  
         
%          if mi_cplmmse(si,i) > mi_cepic_cp(si,i)
%             [mi_cplmmse(si,i) mi_cepic_cp(si,i)]
%             [ser_cplmmse(si,i) ser_cepic_cp(si,i)]
%             [fer_cplmmse(si,i) fer_cepic_cp(si,i)]
%             [mse_lmmsecp_m(si,i) mse_pic_new_1_m(si,i)]
%             []
%          end
         

                        
    end

     Cstamp = clock;
     texttodisp = sprintf('L = %d, SNR = %d, %02d:%02d @ %02d:%02d ended',numPaths,snr(si),Cstamp(2),Cstamp(3),Cstamp(4),Cstamp(5));
     disp(texttodisp)       

    
end

% Pilots overhead
Es = 10*log10(p.M*p.Kon/ ( p.M*p.Kon - (no_observations/nT) )  );


if ~HPC_flag


figure;
semilogy(snr+Es, mean(mse_lmmse.'), '-sm','LineWidth',1)
hold on;
semilogy(snr+Es, mean(mse_lmmsecp.'), '-sb','LineWidth',1)
semilogy(snr+Es, mean(mse_lmmsecp_ofdm.'), '-or','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_1.'), '--o','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_genieD.'), '--ok','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_genieHD.'), '-.ok','LineWidth',1)
semilogy(snr+Es, mean(mse_lmmsez.'), '-db','LineWidth',1) 
semilogy(snr+Es, mean(mse_lmmsecpz.'), '-dm','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1z.'), '-dk','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1z_genie.'), '--dk','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1zcp.'), '-^','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1zcp_genie.'), '--^','LineWidth',1) 
%semilogy(snr, mean(mse_pic_new_1zcp_genie.'), '--om','LineWidth',1) 
title('First Channel Realization')
legend('LMMSE', 'CP-LMMSE','CP-OFDM','CP-LMMSE-PIC','CP-LMMSE-PIC Genie D','CP-LMMSE-PIC Genie H&D','LMMSE Z','CP-LMMSE-Z','LMMSE-Z-PIC','LMMSE-Z-PIC Genie H&D','LMMSE-ZCP-PIC','LMMSE-ZCP-PIC Genie H&D')
grid


    
figure;
semilogy(snr+Es, mean(mse_lmmse_m.'), '-sm','LineWidth',1)
hold on;
semilogy(snr+Es, mean(mse_lmmsecp_m.'), '-sb','LineWidth',1)
semilogy(snr+Es, mean(mse_lmmsecp_ofdm_m.'), '-or','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_1_m.'), '--o','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_genieD_m.'), '--ok','LineWidth',1)
semilogy(snr+Es, mean(mse_pic_new_genieHD_m.'), '-.ok','LineWidth',1)
semilogy(snr+Es, mean(mse_lmmsez_m.'), '-db','LineWidth',1) 
semilogy(snr+Es, mean(mse_lmmsecpz_m.'), '-dm','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1z_m.'), '-dk','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1z_genie_m.'), '--dk','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1zcp_m.'), '-^','LineWidth',1) 
semilogy(snr+Es, mean(mse_pic_new_1zcp_genie_m.'), '--^','LineWidth',1) 
%semilogy(snr, mean(mse_pic_new_1zcp_genie_m.'), '--om','LineWidth',1) 
title('Mean of the time-variant channels')
legend('LMMSE', 'CP-LMMSE','CP-OFDM','CP-LMMSE-PIC','CP-LMMSE-PIC Genie D','CP-LMMSE-PIC Genie H & D','LMMSE Z','CP-LMMSE-Z','LMMSE-Z-PIC','LMMSE-Z-PIC Genie H&D','LMMSE-ZCP-PIC','LMMSE-ZCP-PIC Genie H&D')
grid

if detection_flag
EbNo = snr+Es - 10*log10(p.mu*p.coderate);

figure; 
semilogy(EbNo, mean(bitErr_ofdm.'), '-or','LineWidth',1)
hold on
semilogy(EbNo, mean(bitErr_lmmse.'), '-sm','LineWidth',1)
semilogy(EbNo, mean(bitErr_lmmsez.'), '-db','LineWidth',1)
semilogy(EbNo, mean(bitErr_cplmmse.'), '-sb','LineWidth',1)
semilogy(EbNo, mean(bitErr_cplmmsez.'), '-dm','LineWidth',1)
semilogy(EbNo, mean(bitErr_cepic_cp.'), '-.o','LineWidth',1)
semilogy(EbNo, mean(bitErr_cepicz.'), '-dk','LineWidth',1) 
semilogy(EbNo, mean(bitErr_cepicz_cp.'), '--^k','LineWidth',1) 
grid      
legend('OFDM','LMMSE','LMMSE-Z' , 'CP-LMMSE', 'CP-LMMSE-Z', 'CP-LMMSE-PIC', 'LMMSE-Z-PIC', 'CP-LMMSE-Z')
title('Coded')

figure; 
semilogy(EbNo, mean(fer_ofdm.'), '-or','LineWidth',1)
hold on
semilogy(EbNo, mean(fer_lmmse.'), '-sb','LineWidth',1)
semilogy(EbNo, mean(fer_lmmsez.'), '-db','LineWidth',1)
semilogy(EbNo, mean(fer_cplmmse.'), '-sm','LineWidth',1)
semilogy(EbNo, mean(fer_cplmmsez.'), '-dm','LineWidth',1)
semilogy(EbNo, mean(fer_cepic_cp.'), '-o','LineWidth',1)
semilogy(EbNo, mean(fer_cepicz.'), '-dk','LineWidth',1) 
semilogy(EbNo, mean(fer_cepicz_cp.'), '-^','LineWidth',1) 
grid      
legend('OFDM','LMMSE','SLMMSE' , 'CP-LMMSE', 'CP-SLMMSE', 'CP-LMMSE-PIC', 'SLMMSE-PIC', 'CP-SLMMSE-PIC')
title('Coded FER')

EbNo_uncoded = snr+Es - 10*log10(p.mu);


figure; 
plot(EbNo_uncoded, mean(mi_ofdm.'), '-or','LineWidth',1)
hold on
plot(EbNo_uncoded, mean(mi_lmmse.'), '-sb','LineWidth',1)
plot(EbNo_uncoded, mean(mi_lmmsez.'), '-db','LineWidth',1)
plot(EbNo_uncoded, mean(mi_cplmmse.'), '-sm','LineWidth',1)
plot(EbNo_uncoded, mean(mi_cplmmsez.'), '-dm','LineWidth',1)
plot(EbNo_uncoded, mean(mi_cepic_cp.'), '-o','LineWidth',1)
plot(EbNo_uncoded, mean(mi_cepicz.'), '-dk','LineWidth',1) 
plot(EbNo_uncoded, mean(mi_cepicz_cp.'), '-^','LineWidth',1) 
grid      
legend('OFDM','LMMSE','SLMMSE' , 'CP-LMMSE', 'CP-SLMMSE', 'CP-LMMSE-PIC', 'SLMMSE-PIC', 'CP-SLMMSE-PIC')
title('Mutual Information')

figure; 
semilogy(EbNo_uncoded, mean(ser_ofdm.'), '-or','LineWidth',1)
hold on
semilogy(EbNo_uncoded, mean(ser_lmmse.'), '-sm','LineWidth',1)
semilogy(EbNo_uncoded, mean(ser_lmmsez.'), '-db','LineWidth',1)
semilogy(EbNo_uncoded, mean(ser_cplmmse.'), '-sb','LineWidth',1)
semilogy(EbNo_uncoded, mean(ser_cplmmsez.'), '-dm','LineWidth',1)
semilogy(EbNo_uncoded, mean(ser_cepic_cp.'), '-.o','LineWidth',1)
semilogy(EbNo_uncoded, mean(ser_cepicz.'), '-dk','LineWidth',1) 
semilogy(EbNo_uncoded, mean(ser_cepicz_cp.'), '--^k','LineWidth',1) 
grid      
legend('OFDM','LMMSE','LMMSE-Z' , 'CP-LMMSE', 'CP-LMMSE-Z', 'CP-LMMSE-PIC', 'LMMSE-Z-PIC', 'CP-LMMSE-Z')
title('Uncoded')
end

else
    if fd < 1
    savestring = sprintf('a38_CP_aided_paper_seed%d_MCS%d.mat',rnd_seed,MCS);
    %Folder = 'C:\Documents\MATLAB\GFDM\savedWorkspaces';
    save(['./a38_MCS3_fd70/' savestring])
    else
    savestring = sprintf('a38_CP_aided_paper_seed%d_MCS%d_fd%d.mat',rnd_seed,MCS,fd);
    %Folder = 'C:\Documents\MATLAB\GFDM\savedWorkspaces';
    save(['./a38_MCS3_fd70/' savestring])
    end
    exit
end



