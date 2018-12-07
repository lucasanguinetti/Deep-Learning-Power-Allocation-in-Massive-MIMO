clear;

Input_tr = [];

%Prepare MR maxprod
Output_tr_MR_maxmin = [];
Output_tr_MR_maxmin_cell_1 = [];
Output_tr_MR_maxmin_cell_2 = [];
Output_tr_MR_maxmin_cell_3 = [];
Output_tr_MR_maxmin_cell_4 = [];

%Prepare MR maxprod
Output_tr_MMMSE_maxmin = [];
Output_tr_MMMSE_maxmin_cell_1 = [];
Output_tr_MMMSE_maxmin_cell_2 = [];
Output_tr_MMMSE_maxmin_cell_3 = [];
Output_tr_MMMSE_maxmin_cell_4 = [];


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
    output_MR_maxmin_all_cells = reshape(output_MR_maxmin(:,:,:), K*L,[]);
    %     output_MR_maxmin_all_cells = [output_MR_maxmin_all_cells; sum(output_MR_maxmin_all_cells,1)];
    
    output_MR_maxmin_cell_1 = reshape(output_MR_maxmin(:,1,:), K,[]);
    %     output_MR_maxmin_cell_1 = [output_MR_maxmin_cell_1; sum(output_MR_maxmin_cell_1,1)];
    
    output_MR_maxmin_cell_2 = reshape(output_MR_maxmin(:,2,:), K,[]);
    %     output_MR_maxmin_cell_2 = [output_MR_maxmin_cell_2; sum(output_MR_maxmin_cell_2,1)];
    
    output_MR_maxmin_cell_3 = reshape(output_MR_maxmin(:,3,:), K,[]);
    %     output_MR_maxmin_cell_3 = [output_MR_maxmin_cell_3; sum(output_MR_maxmin_cell_3,1)];
    
    output_MR_maxmin_cell_4 = reshape(output_MR_maxmin(:,4,:), K,[]);
    %     output_MR_maxmin_cell_4 = [output_MR_maxmin_cell_4; sum(output_MR_maxmin_cell_4,1)];
    
    % Output for NN with MMMSE
    output_MMMSE_maxmin_all_cells = reshape(output_MMMSE_maxmin(:,:,:), K*L,[]);
    
    output_MMMSE_maxmin_cell_1 = reshape(output_MMMSE_maxmin(:,1,:), K,[]);
    %     output_MMMSE_maxmin_cell_1 = [output_MMMSE_maxmin_cell_1; sum(output_MMMSE_maxmin_cell_1,1)];
    
    output_MMMSE_maxmin_cell_2 = reshape(output_MMMSE_maxmin(:,2,:), K,[]);
    %     output_MMMSE_maxmin_cell_2 = [output_MMMSE_maxmin_cell_2; sum(output_MMMSE_maxmin_cell_2,1)];
    
    output_MMMSE_maxmin_cell_3 = reshape(output_MMMSE_maxmin(:,3,:), K,[]);
    %     output_MMMSE_maxmin_cell_3 = [output_MMMSE_maxmin_cell_3; sum(output_MMMSE_maxmin_cell_3,1)];
    
    output_MMMSE_maxmin_cell_4 = reshape(output_MMMSE_maxmin(:,4,:), K,[]);
    %     output_MMMSE_maxmin_cell_4 = [output_MMMSE_maxmin_cell_4; sum(output_MMMSE_maxmin_cell_4,1)];
    
    % Concatenate input data for NN
    Input_tr = [Input_tr, input]; %#ok<*AGROW>
    
    % Concatenate output data for NN
    Output_tr_MR_maxmin = [Output_tr_MR_maxmin, output_MR_maxmin_all_cells];
    
    Output_tr_MR_maxmin_cell_1 = [Output_tr_MR_maxmin_cell_1, output_MR_maxmin_cell_1];
    Output_tr_MR_maxmin_cell_2 = [Output_tr_MR_maxmin_cell_2, output_MR_maxmin_cell_2];
    Output_tr_MR_maxmin_cell_3 = [Output_tr_MR_maxmin_cell_3, output_MR_maxmin_cell_3];
    Output_tr_MR_maxmin_cell_4 = [Output_tr_MR_maxmin_cell_4, output_MR_maxmin_cell_4];
    
    
    % Concatenate output data for NN with M-MMSE
    Output_tr_MMMSE_maxmin = [Output_tr_MMMSE_maxmin, output_MMMSE_maxmin_all_cells];
    
    Output_tr_MMMSE_maxmin_cell_1 = [Output_tr_MMMSE_maxmin_cell_1, output_MMMSE_maxmin_cell_1];
    Output_tr_MMMSE_maxmin_cell_2 = [Output_tr_MMMSE_maxmin_cell_2, output_MMMSE_maxmin_cell_2];
    Output_tr_MMMSE_maxmin_cell_3 = [Output_tr_MMMSE_maxmin_cell_3, output_MMMSE_maxmin_cell_3];
    Output_tr_MMMSE_maxmin_cell_4 = [Output_tr_MMMSE_maxmin_cell_4, output_MMMSE_maxmin_cell_4];
    
end


% Conversion in dB for input data
Input_tr_dB = 10*log10(Input_tr);

Input_tr_normalized = (Input_tr - mean(Input_tr,2))./std(Input_tr,0,2);

Input_tr_dB_normalized = (Input_tr_dB - mean(Input_tr_dB,2))./std(Input_tr_dB,0,2);


Output_tr_MR_maxmin_normalized_cell_1 = (Output_tr_MR_maxmin_cell_1 - mean(Output_tr_MR_maxmin_cell_1,2))./std(Output_tr_MR_maxmin_cell_1,0,2);

Output_tr_MR_maxmin_dB_cell_1 = 10*log(Output_tr_MR_maxmin_cell_1);

Output_tr_MR_maxmin_dB_normalized_cell_1 = (Output_tr_MR_maxmin_dB_cell_1 - mean(Output_tr_MR_maxmin_dB_cell_1,2))./std(Output_tr_MR_maxmin_dB_cell_1,0,2);


Output_tr_MR_maxmin_normalized_cell_2 = (Output_tr_MR_maxmin_cell_2 - mean(Output_tr_MR_maxmin_cell_2,2))./std(Output_tr_MR_maxmin_cell_2,0,2);

Output_tr_MR_maxmin_dB_cell_2 = 10*log(Output_tr_MR_maxmin_cell_2);

Output_tr_MR_maxmin_dB_normalized_cell_2 = (Output_tr_MR_maxmin_dB_cell_2 - mean(Output_tr_MR_maxmin_dB_cell_2,2))./std(Output_tr_MR_maxmin_dB_cell_2,0,2);


Output_tr_MR_maxmin_normalized_cell_3 = (Output_tr_MR_maxmin_cell_3 - mean(Output_tr_MR_maxmin_cell_3,2))./std(Output_tr_MR_maxmin_cell_3,0,2);

Output_tr_MR_maxmin_dB_cell_3 = 10*log(Output_tr_MR_maxmin_cell_3);

Output_tr_MR_maxmin_dB_normalized_cell_3 = (Output_tr_MR_maxmin_dB_cell_3 - mean(Output_tr_MR_maxmin_dB_cell_3,2))./std(Output_tr_MR_maxmin_dB_cell_3,0,2);


Output_tr_MR_maxmin_normalized_cell_4 = (Output_tr_MR_maxmin_cell_4 - mean(Output_tr_MR_maxmin_cell_4,2))./std(Output_tr_MR_maxmin_cell_4,0,2);

Output_tr_MR_maxmin_dB_cell_4 = 10*log(Output_tr_MR_maxmin_cell_4);

Output_tr_MR_maxmin_dB_normalized_cell_4 = (Output_tr_MR_maxmin_dB_cell_4 - mean(Output_tr_MR_maxmin_dB_cell_4,2))./std(Output_tr_MR_maxmin_dB_cell_4,0,2);


Output_tr_MMMSE_maxmin_normalized_cell_1 = (Output_tr_MMMSE_maxmin_cell_1 - mean(Output_tr_MMMSE_maxmin_cell_1,2))./std(Output_tr_MMMSE_maxmin_cell_1,0,2);

Output_tr_MMMSE_maxmin_dB_cell_1 = 10*log(Output_tr_MMMSE_maxmin_cell_1);

Output_tr_MMMSE_maxmin_dB_normalized_cell_1 = (Output_tr_MMMSE_maxmin_dB_cell_1 - mean(Output_tr_MMMSE_maxmin_dB_cell_1,2))./std(Output_tr_MMMSE_maxmin_dB_cell_1,0,2);


Output_tr_MMMSE_maxmin_normalized_cell_2 = (Output_tr_MMMSE_maxmin_cell_2 - mean(Output_tr_MMMSE_maxmin_cell_2,2))./std(Output_tr_MMMSE_maxmin_cell_2,0,2);

Output_tr_MMMSE_maxmin_dB_cell_2 = 10*log(Output_tr_MMMSE_maxmin_cell_2);

Output_tr_MMMSE_maxmin_dB_normalized_cell_2 = (Output_tr_MMMSE_maxmin_dB_cell_2 - mean(Output_tr_MMMSE_maxmin_dB_cell_2,2))./std(Output_tr_MMMSE_maxmin_dB_cell_2,0,2);


Output_tr_MMMSE_maxmin_normalized_cell_3 = (Output_tr_MMMSE_maxmin_cell_3 - mean(Output_tr_MMMSE_maxmin_cell_3,2))./std(Output_tr_MMMSE_maxmin_cell_3,0,2);

Output_tr_MMMSE_maxmin_dB_cell_3 = 10*log(Output_tr_MMMSE_maxmin_cell_3);

Output_tr_MMMSE_maxmin_dB_normalized_cell_3 = (Output_tr_MMMSE_maxmin_dB_cell_3 - mean(Output_tr_MMMSE_maxmin_dB_cell_3,2))./std(Output_tr_MMMSE_maxmin_dB_cell_3,0,2);


Output_tr_MMMSE_maxmin_normalized_cell_4 = (Output_tr_MMMSE_maxmin_cell_4 - mean(Output_tr_MMMSE_maxmin_cell_4,2))./std(Output_tr_MMMSE_maxmin_cell_4,0,2);

Output_tr_MMMSE_maxmin_dB_cell_4 = 10*log(Output_tr_MMMSE_maxmin_cell_4);

Output_tr_MMMSE_maxmin_dB_normalized_cell_4 = (Output_tr_MMMSE_maxmin_dB_cell_4 - mean(Output_tr_MMMSE_maxmin_dB_cell_4,2))./std(Output_tr_MMMSE_maxmin_dB_cell_4,0,2);


figure;
hold on; box on;
histogram(Output_tr_MR_maxmin_dB_cell_1(:))
histogram(Output_tr_MR_maxmin_dB_normalized_cell_1(:))


save('dataset_maxmin.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized',...
    'Output_tr_MR_maxmin_cell_1','Output_tr_MR_maxmin_dB_cell_2',...
    'Output_tr_MR_maxmin_cell_3','Output_tr_MR_maxmin_dB_cell_4',...
    'Output_tr_MR_maxmin_normalized_cell_1','Output_tr_MR_maxmin_normalized_cell_2',...
    'Output_tr_MR_maxmin_normalized_cell_3','Output_tr_MR_maxmin_normalized_cell_4',...
    'Output_tr_MR_maxmin_dB_normalized_cell_1','Output_tr_MR_maxmin_dB_normalized_cell_2',...
    'Output_tr_MR_maxmin_dB_normalized_cell_3','Output_tr_MR_maxmin_dB_normalized_cell_4',...
    'Output_tr_MMMSE_maxmin_cell_1','Output_tr_MMMSE_maxmin_cell_2',...
    'Output_tr_MMMSE_maxmin_cell_3','Output_tr_MMMSE_maxmin_cell_4',...
    'Output_tr_MMMSE_maxmin_normalized_cell_1','Output_tr_MMMSE_maxmin_normalized_cell_2',...
    'Output_tr_MMMSE_maxmin_normalized_cell_3','Output_tr_MMMSE_maxmin_normalized_cell_4',...
    'Output_tr_MMMSE_maxmin_dB_normalized_cell_1','Output_tr_MMMSE_maxmin_dB_normalized_cell_2',...
    'Output_tr_MMMSE_maxmin_dB_normalized_cell_3','Output_tr_MMMSE_maxmin_dB_normalized_cell_4');
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

save('testset_maxmin.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized');

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


save('testset_maxmin.mat','Pmax', 'Input_tr','Input_tr_dB','Input_tr_normalized','Input_tr_dB_normalized',...
    'Input_tr_cell1','Input_tr_dB_cell1','Input_tr_dB_normalized_cell1',...
    'Input_tr_cell2','Input_tr_dB_cell2','Input_tr_dB_normalized_cell2',...
    'Input_tr_cell3','Input_tr_dB_cell3','Input_tr_dB_normalized_cell3',...
    'Input_tr_cell4','Input_tr_dB_cell4','Input_tr_dB_normalized_cell4');