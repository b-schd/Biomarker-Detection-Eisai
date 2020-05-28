function params = initialize_task

% Specify toolbox paths
% addpath(genpath('../ieeg-matlab-1.13.2'));
% addpath(genpath('../portal-matlab-tools/Analysis'))
% addpath(genpath('../portal-matlab-tools/Utilities'))


% Initialize variables
%params.dataKey = humanNV_dataKey
params.annotFile='/Users/bscheid/Documents/LittLab/PROJECTS/p07_PIG-PTE/test_annots.xlsx';

params.datasetID = {
    'PIG_PTE',
    };

%% includes emg channels
% params.groupChannels = {
%    {'EEG EEG1.1A-B','EEG EEG2.1A-B','EMG EMG.1'},  
%    {'EEG EEG1.2A-B','EEG EEG2.2A-B','EMG EMG.2'},
%    %{'EEG EEG1.3A-B','EEG EEG2.3A-B','EMG EMG.3'}, %not used
%    {'EEG EEG1A-B','EEG EEG2A-B','EMG EMG'},
%    };

%%
params.groupChannels = {
    {'EEG EEG1.1A-B','EEG EEG2.1A-B'},  
    {'EEG EEG1.2A-B','EEG EEG2.2A-B'},
    {'EEG EEG1.3A-B','EEG EEG2.3A-B'}, %not used
    {'EEG EEG1A-B','EEG EEG2A-B'},
    };


params.IEEGid = 'bscheid';
params.IEEGpwd = '/Users/bscheid/ieeg_pwd.bin';
params.spike_layer = '';
params.marked_seizure_layers = {'True_Seizures', 'Non_Seizures'};
