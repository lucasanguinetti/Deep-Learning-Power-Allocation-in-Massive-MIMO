clear;

Input_tr = [];

%Prepare MR maxprod
Output_tr_MR_maxprod = [];
Output_tr_MR_maxprod_cell_1 = [];
Output_tr_MR_maxprod_cell_2 = [];
Output_tr_MR_maxprod_cell_3 = [];
Output_tr_MR_maxprod_cell_4 = [];

%Prepare MR maxprod
Output_tr_MMMSE_maxprod = [];
Output_tr_MMMSE_maxprod_cell_1 = [];
Output_tr_MMMSE_maxprod_cell_2 = [];
Output_tr_MMMSE_maxprod_cell_3 = [];
Output_tr_MMMSE_maxprod_cell_4 = [];


% BS positions
BS_positions = [250+1i*250; 250+1i*750; 750+1i*250; 750+1i*750];

indexes = 1:2; 

for kk = 1:length(indexes)
    
    ii = indexes(kk);
    
    % Load data file
    load(['MyDataFile_' mat2str(ii) '.mat'])
    
    % Position of UEs in the network
    input_positions_reshaped = reshape(input_positions,K*L,[]);
    
    % Prepare BS positions
    BS_positions_rep = repmat(BS_positions,1,size(input_positions_reshaped,2));
    
    % input for NN
%     input = [real(input_positions_reshaped); imag(input_positions_reshaped); real(BS_positions_rep); imag(BS_positions_rep)];
    input = [real(input_positions_reshaped); imag(input_positions_reshaped)];
    
    % Output for NN with MR and MMMSE for all cells
    output_MR_maxprod_all_cells = reshape(output_MR_maxprod(:,:,:), K*L,[]);
    output_MR_maxprod_all_cells = [output_MR_maxprod_all_cells; sum(output_MR_maxprod_all_cells,1)];
    
    output_MR_maxprod_cell_1 = reshape(output_MR_maxprod(:,1,:), K,[]);
    output_MR_maxprod_cell_1 = [output_MR_maxprod_cell_1; sum(output_MR_maxprod_cell_1,1)];
    
    output_MR_maxprod_cell_2 = reshape(output_MR_maxprod(:,2,:), K,[]);
    output_MR_maxprod_cell_2 = [output_MR_maxprod_cell_2; sum(output_MR_maxprod_cell_2,1)];
    
    output_MR_maxprod_cell_3 = reshape(output_MR_maxprod(:,3,:), K,[]);
    output_MR_maxprod_cell_3 = [output_MR_maxprod_cell_3; sum(output_MR_maxprod_cell_3,1)];
    
    output_MR_maxprod_cell_4 = reshape(output_MR_maxprod(:,4,:), K,[]);
    output_MR_maxprod_cell_4 = [output_MR_maxprod_cell_4; sum(output_MR_maxprod_cell_4,1)];
    
    % Output for NN with MMMSE
    output_MMMSE_maxprod_all_cells = reshape(output_MMMSE_maxprod(:,:,:), K*L,[]);
    
    output_MMMSE_maxprod_cell_1 = reshape(output_MMMSE_maxprod(:,1,:), K,[]);
    output_MMMSE_maxprod_cell_1 = [output_MMMSE_maxprod_cell_1; sum(output_MMMSE_maxprod_cell_1,1)];
    
    output_MMMSE_maxprod_cell_2 = reshape(output_MMMSE_maxprod(:,2,:), K,[]);
    output_MMMSE_maxprod_cell_2 = [output_MMMSE_maxprod_cell_2; sum(output_MMMSE_maxprod_cell_2,1)];
    
    output_MMMSE_maxprod_cell_3 = reshape(output_MMMSE_maxprod(:,3,:), K,[]);
    output_MMMSE_maxprod_cell_3 = [output_MMMSE_maxprod_cell_3; sum(output_MMMSE_maxprod_cell_3,1)];
    
    output_MMMSE_maxprod_cell_4 = reshape(output_MMMSE_maxprod(:,4,:), K,[]);
    output_MMMSE_maxprod_cell_4 = [output_MMMSE_maxprod_cell_4; sum(output_MMMSE_maxprod_cell_4,1)];
    
    % Concatenate input data for NN
    Input_tr = [Input_tr, input]; %#ok<*AGROW>
    
    % Concatenate output data for NN
    Output_tr_MR_maxprod = [Output_tr_MR_maxprod, output_MR_maxprod_all_cells];
    
    Output_tr_MR_maxprod_cell_1 = [Output_tr_MR_maxprod_cell_1, output_MR_maxprod_cell_1];
    Output_tr_MR_maxprod_cell_2 = [Output_tr_MR_maxprod_cell_2, output_MR_maxprod_cell_2];
    Output_tr_MR_maxprod_cell_3 = [Output_tr_MR_maxprod_cell_3, output_MR_maxprod_cell_3];
    Output_tr_MR_maxprod_cell_4 = [Output_tr_MR_maxprod_cell_4, output_MR_maxprod_cell_4];
    
    
    % Concatenate output data for NN with M-MMSE
    Output_tr_MMMSE_maxprod = [Output_tr_MMMSE_maxprod, output_MMMSE_maxprod_all_cells];
    
    Output_tr_MMMSE_maxprod_cell_1 = [Output_tr_MMMSE_maxprod_cell_1, output_MMMSE_maxprod_cell_1];
    Output_tr_MMMSE_maxprod_cell_2 = [Output_tr_MMMSE_maxprod_cell_2, output_MMMSE_maxprod_cell_2];
    Output_tr_MMMSE_maxprod_cell_3 = [Output_tr_MMMSE_maxprod_cell_3, output_MMMSE_maxprod_cell_3];
    Output_tr_MMMSE_maxprod_cell_4 = [Output_tr_MMMSE_maxprod_cell_4, output_MMMSE_maxprod_cell_4];
    
end

% Conversion in dB for input data
Input_tr_dB = 10*log10(Input_tr);

Input_tr_normalized = (Input_tr - mean(Input_tr,2))./std(Input_tr,0,2);

Input_tr_dB_normalized = (Input_tr_dB - mean(Input_tr_dB,2))./std(Input_tr_dB,0,2);

% % Normalized
% Input_tr = 10*log10(Input_tr);
% Input_tr_normalized = zeros(size(Input_tr));
% 
% Input_tr_normalized(1:end-2*L,:) = (Input_tr(1:end-2*L,:) - mean(Input_tr(1:end-2*L,:),2))./std(Input_tr(1:end-2*L,:),0,2);
% 
% Input_tr_normalized(end-2*L+1:end,:) = Input_tr(end-2*L+1:end,:)./max(Input_tr(end-2*L+1:end,:),[],2);


save('dataset_maxprod.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized',...
    'Output_tr_MR_maxprod',...
    'Output_tr_MR_maxprod_cell_1','Output_tr_MR_maxprod_cell_2',...
    'Output_tr_MR_maxprod_cell_3','Output_tr_MR_maxprod_cell_4',...
    'Output_tr_MMMSE_maxprod',...
    'Output_tr_MMMSE_maxprod_cell_1','Output_tr_MMMSE_maxprod_cell_2',...
    'Output_tr_MMMSE_maxprod_cell_3','Output_tr_MMMSE_maxprod_cell_4');

clear;

load('MyDataFile_0.mat')

% BS positions
BS_positions = [250+1i*250; 250+1i*750; 750+1i*250; 750+1i*750];

% Position of UEs in the network
input_positions_reshaped = reshape(input_positions,K*L,[]);
% Prepare BS positions
BS_positions_rep = repmat(BS_positions,1,size(input_positions_reshaped,2));
% input for NN
Input_tr = [real(input_positions_reshaped); imag(input_positions_reshaped)];%; real(BS_positions_rep); imag(BS_positions_rep)];

% Conversion in dB for input data
Input_tr_dB = 10*log10(Input_tr);

Input_tr_normalized = (Input_tr - mean(Input_tr,2))./std(Input_tr,0,2);

Input_tr_dB_normalized = (Input_tr_dB - mean(Input_tr_dB,2))./std(Input_tr_dB,0,2);

save('testset_maxprod.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized');

%
input_positions_cell1 = squeeze(input_positions(:,1,:));
input_positions_cell2 = squeeze(input_positions(:,2,:));
input_positions_cell3 = squeeze(input_positions(:,3,:));
input_positions_cell4 = squeeze(input_positions(:,4,:));

% Prepare BS positions

%Cell 1
input_positions_reshaped = [input_positions_cell1;input_positions_cell2;input_positions_cell3;input_positions_cell4];

Input_tr_cell1 = [real(input_positions_reshaped); imag(input_positions_reshaped)];

Input_tr_normalized_cell1 = (Input_tr_cell1 - mean(Input_tr_cell1,2))./std(Input_tr_cell1,0,2);

Input_tr_dB_cell1 = 10*log10(Input_tr_cell1);

Input_tr_dB_normalized_cell1 = (Input_tr_dB_cell1 - mean(Input_tr_dB_cell1,2))./std(Input_tr_dB_cell1,0,2);

%Cell 2
input_positions_reshaped = [input_positions_cell2;input_positions_cell3;input_positions_cell4;input_positions_cell1];

Input_tr_cell2 = [real(input_positions_reshaped); imag(input_positions_reshaped)];

Input_tr_normalized_cell3 = (Input_tr_cell2 - mean(Input_tr_cell2,2))./std(Input_tr_cell2,0,2);

Input_tr_dB_cell2 = 10*log10(Input_tr_cell2);

Input_tr_dB_normalized_cell2 = (Input_tr_dB_cell2 - mean(Input_tr_dB_cell2,2))./std(Input_tr_dB_cell2,0,2);

%Cell 3
input_positions_reshaped = [input_positions_cell3;input_positions_cell4;input_positions_cell1;input_positions_cell2];

Input_tr_cell3 = [real(input_positions_reshaped); imag(input_positions_reshaped)];

Input_tr_normalized_cell3 = (Input_tr_cell3 - mean(Input_tr_cell3,2))./std(Input_tr_cell3,0,2);

Input_tr_dB_cell3 = 10*log10(Input_tr_cell3);

Input_tr_dB_normalized_cell3 = (Input_tr_dB_cell3 - mean(Input_tr_dB_cell3,2))./std(Input_tr_dB_cell3,0,2);

%Cell 4
input_positions_reshaped = [input_positions_cell4;input_positions_cell1;input_positions_cell2;input_positions_cell3];

Input_tr_cell4 = [real(input_positions_reshaped); imag(input_positions_reshaped)];

Input_tr_normalized_cell4 = (Input_tr_cell4 - mean(Input_tr_cell4,2))./std(Input_tr_cell4,0,2);

Input_tr_dB_cell4 = 10*log10(Input_tr_cell4);

Input_tr_dB_normalized_cell4 = (Input_tr_dB_cell4 - mean(Input_tr_dB_cell4,2))./std(Input_tr_dB_cell4,0,2);


save('testset_maxprod.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized',...
    'Input_tr_cell1','Input_tr_dB_cell1','Input_tr_dB_normalized_cell1',...
    'Input_tr_cell2','Input_tr_dB_cell2','Input_tr_dB_normalized_cell2',...
    'Input_tr_cell3','Input_tr_dB_cell3','Input_tr_dB_normalized_cell3',...
    'Input_tr_cell4','Input_tr_dB_cell4','Input_tr_dB_normalized_cell4');