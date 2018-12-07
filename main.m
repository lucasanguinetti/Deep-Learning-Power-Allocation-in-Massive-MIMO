%This Matlab script can be used to reproduce Figure 7.2 in the monograph:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.0 (Last edited: 2017-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.
%
%Note: This script require additional software packages to be used, which
%need to be downloaded and installed separately. These packages are
%developed independently and are delivered with separate licenses.
%
%The downlink power allocation is optimized using CVX from CVX Research,
%Inc. (http://cvxr.com/cvx/). This script has been tested with CVX 2.1,
%using the solver Mosek, version 7.1.0.12. We discourage the use of the
%solvers SDPT3 and SeDuMi since these crashed during the test.


%Empty workspace and close figures
close all;
clear;
       
rng('shuffle')

%Number of BSs
L = 4;

%Number of UEs per BS
K = 5;

%Number of BS antennas
M = 100;

%Define the pilot reuse factor
f = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 5000;

%Select the number of channel realizations per setup
nbrOfRealizations = 100;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per UE (mW)
rho = 100;

%Maximum downlink transmit power per BS (mW)
Pmax = K*rho;

%Compute downlink power per UE in case of equal power allocation
rhoEqual = (Pmax/K)*ones(K,L);

%Define noise figure at BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
SE_MR_equal = zeros(K,L,nbrOfSetups);
SE_RZF_equal = zeros(K,L,nbrOfSetups);
SE_MMMSE_equal = zeros(K,L,nbrOfSetups);

SE_MR_maxmin = zeros(K,L,nbrOfSetups);
SE_RZF_maxmin = zeros(K,L,nbrOfSetups);
SE_MMMSE_maxmin = zeros(K,L,nbrOfSetups);

SE_MR_maxprod = zeros(K,L,nbrOfSetups);
SE_RZF_maxprod = zeros(K,L,nbrOfSetups);
SE_MMMSE_maxprod = zeros(K,L,nbrOfSetups);


% Prepare to save data input for ANN with MR, RZF, MMSE
input_channelGaindB = zeros(K,L,L,nbrOfSetups);
input_positions = zeros(K,L,nbrOfSetups);

output_MR_maxmin = zeros(K,L,nbrOfSetups);
output_RZF_maxmin = zeros(K,L,nbrOfSetups);
output_MMMSE_maxmin = zeros(K,L,nbrOfSetups);

output_MR_maxprod = zeros(K,L,nbrOfSetups);
output_RZF_maxprod = zeros(K,L,nbrOfSetups);
output_MMMSE_maxprod = zeros(K,L,nbrOfSetups);

% Prepare to save average signal power for ANN with MR, RZF, MMSE
signal_MR_AAN = zeros(K,L,nbrOfSetups);
signal_RZF_AAN = zeros(K,L,nbrOfSetups);
signal_MMMSE_AAN = zeros(K,L,nbrOfSetups);
% Prepare to save average interference power for ANN with MR, RZF, MMSE
interf_MR_AAN = zeros(K,L,K,L,nbrOfSetups);
interf_RZF_AAN = zeros(K,L,K,L,nbrOfSetups);
interf_MMMSE_AAN = zeros(K,L,K,L,nbrOfSetups);




%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB,UEpositions] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    % Saving input data for ANN with all schemes
    input_channelGaindB(:,:,:,n) = channelGaindB;
    
        % Saving input data for ANN with all schemes
    input_positions(:,:,n) = UEpositions;
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all number of antennas
    for m = 1:length(M)
        
        
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
        
        %Compute the signal and interference terms of the DL SEs using the
        %hardening bound in Theorem 4.6
        [signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE,interf_MMMSE,prelogFactor] = functionComputeSINR_DL(H,Hhat,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
        
        %Delete large matrices
        clear Hhat C Rscaled;
        
        % Saving average signal power for ANN with MR, RZF, MMSE
        signal_MR_AAN(:,:,n) = signal_MR; 
        signal_RZF_AAN(:,:,n) = signal_RZF; 
        signal_MMMSE_AAN(:,:,n) = signal_MMMSE;
        
        
        % Saving interference power for ANN with MR, RZF, MMSE
        interf_MR_AAN(:,:,:,:,n) = interf_MR; 
        interf_RZF_AAN(:,:,:,:,n) = interf_RZF; 
        interf_MMMSE_AAN(:,:,:,:,n) = interf_MMMSE;
        
        
        %Compute the SEs with equal power allocation
        SE_MR_equal(:,:,n) = functionComputeSE_DL_poweralloc(rhoEqual,signal_MR,interf_MR,prelogFactor);
        SE_RZF_equal(:,:,n) = functionComputeSE_DL_poweralloc(rhoEqual,signal_RZF,interf_RZF,prelogFactor);
        SE_MMMSE_equal(:,:,n) = functionComputeSE_DL_poweralloc(rhoEqual,signal_MMMSE,interf_MMMSE,prelogFactor);
        
        %Compute the SEs with max-min fairness power allocation
        disp('Solve max-min fairness problem'); %Output simulation progress
        [SE_MR_maxmin(:,:,n), output_MR_maxmin(:,:,n)] = functionPowerOptimization_maxmin(signal_MR,interf_MR,Pmax,prelogFactor);
        [SE_RZF_maxmin(:,:,n), output_RZF_maxmin(:,:,n)] = functionPowerOptimization_maxmin(signal_RZF,interf_RZF,Pmax,prelogFactor);
        [SE_MMMSE_maxmin(:,:,n), output_MMMSE_maxmin(:,:,n)] = functionPowerOptimization_maxmin(signal_MMMSE,interf_MMMSE,Pmax,prelogFactor);
        
        %Compute the SEs with max-min fairness power allocation
        disp('Solve max product SINR problem'); %Output simulation progress
        [SE_MR_maxprod(:,:,n), output_MR_maxprod(:,:,n)]  = functionPowerOptimization_prodSINR(signal_MR,interf_MR,Pmax,prelogFactor);
        [SE_RZF_maxprod(:,:,n), output_RZF_maxprod(:,:,n)] = functionPowerOptimization_prodSINR(signal_RZF,interf_RZF,Pmax,prelogFactor);
        [SE_MMMSE_maxprod(:,:,n), output_MMMSE_maxprod(:,:,n)] = functionPowerOptimization_prodSINR(signal_MMMSE,interf_MMMSE,Pmax,prelogFactor);
        
        
    end
    
    %Delete large matrices
    clear R;
    
        
    % Save every 1000 setups
    if mod(n,1000) == 0
       
        save('MyDataFile.mat', 'L', 'K', 'M', 'Pmax', 'prelogFactor',...
            'input_channelGaindB', 'input_positions', 'nbrOfSetups', 'nbrOfRealizations',...
            'signal_MR_AAN', 'signal_RZF_AAN', 'signal_MMMSE_AAN',...
            'interf_MR_AAN', 'interf_RZF_AAN', 'interf_MMMSE_AAN',...
            'output_MR_maxmin', 'output_RZF_maxmin', 'output_MMMSE_maxmin',...
            'output_MR_maxprod', 'output_RZF_maxprod', 'output_MMMSE_maxprod',...
            'SE_MR_equal', 'SE_RZF_equal', 'SE_MMMSE_equal',...
            'SE_MR_maxmin', 'SE_RZF_maxmin', 'SE_MMMSE_maxmin',...
            'SE_MR_maxprod', 'SE_RZF_maxprod', 'SE_MMMSE_maxprod');
        
    end
    
end


%% Plot the simulation results
% 
% figure(1);
% hold on; box on;
% 
% plot(sort(SE_MR_maxmin(:)),linspace(0,1,K*L*nbrOfSetups),'k--','LineWidth',1);
% plot(sort(SE_MR_equal(:)),linspace(0,1,K*L*nbrOfSetups),'k-','LineWidth',1);
% plot(sort(SE_MR_maxprod(:)),linspace(0,1,K*L*nbrOfSetups),'k-.','LineWidth',1);
% 
% plot(sort(SE_RZF_equal(:)),linspace(0,1,K*L*nbrOfSetups),'b-','LineWidth',1);
% plot(sort(SE_MMMSE_equal(:)),linspace(0,1,K*L*nbrOfSetups),'r','LineWidth',1);
% 
% plot(sort(SE_RZF_maxprod(:)),linspace(0,1,K*L*nbrOfSetups),'b-.','LineWidth',1);
% plot(sort(SE_MMMSE_maxprod(:)),linspace(0,1,K*L*nbrOfSetups),'r-.','LineWidth',1);
% 
% plot(sort(SE_RZF_maxmin(:)),linspace(0,1,K*L*nbrOfSetups),'b--','LineWidth',1);
% plot(sort(SE_MMMSE_maxmin(:)),linspace(0,1,K*L*nbrOfSetups),'r--','LineWidth',1);
% 
% legend('Max-min fairness','Equal power','Max product SINR','Location','SouthEast');
% 
% xlabel('SE per UE [bit/s/Hz]');
% ylabel('CDF');
