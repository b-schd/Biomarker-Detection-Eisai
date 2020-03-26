function params = initialize_task

% Initialize variables
%params.dataKey = humanNV_dataKey

params.datasetID = {
    'JAW_72',
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
params.IEEGpwd = 'bscheid_ieeglogin.bin';
params.spike_layer = '';
params.marked_seizure_layer = 'True_Seizures';
