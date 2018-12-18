function params = initialize_task

% Initialize variables
%params.dataKey = humanNV_dataKey

params.datasetID = {
    'Jensen_Eisai_groupA',
   'Jensen_Eisai_groupB',
   'Jensen_Eisai_groupC',
  'Jensen_Eisai_groupD',
   'Jensen_Eisai_groupE',
   'Jensen_Eisai_groupF',
   'Jensen_Eisai_groupG',
   'Jensen_Eisai_groupH',
   'Jensen_Eisai_groupI',
   'Jensen_Eisai_groupJ',
    };

params.groupChannels = {
    {'EEG EEG1.1A-B','EEG EEG2.1A-B','EMG EMG.1'},  
    {'EEG EEG1.2A-B','EEG EEG2.2A-B','EMG EMG.2'},
    %{'EEG EEG1.3A-B','EEG EEG2.3A-B','EMG EMG.3'}, %not used
    {'EEG EEG1A-B','EEG EEG2A-B','EMG EMG'},
    };


params.IEEGid = 'hoameng';
params.IEEGpwd = 'hoa_ieeglogin.bin';
params.spike_layer = '';
params.marked_seizure_layer = 'True_Seizures';
%Input:
%   'dataset'       :   IEEG dataset
%   'channels'      :   vector of channels
%   'winLen'        :   winLen = vector of windowlengths (s) to calculate
%                       features over
%   'outlabel'      :   Suffix to save features to
%   (datasetname_outlabel.mat)
%   'filtFlag'      :   [0/1] 1: filter data [1 70] bandpass, [58 62]
%   bandstop
%   'filtCheck'     :   [0/1] 1: Plot data before and after filtering for
%   one block for manual checking
%Output:
%   'feat'          :   Calculated feature for each window
params.feature = 'll';
params.winLen = 60;%s
params.winDisp = 60;
params.filtFlag = 1;
params.diag = 0;
params.blockLen = 1*60*60;
params.channels = 1:16;
params.saveLabel = 'LL';
params.timeOfInterest = [0 60*60*24*100]; %[start stop](s)
% 
% switch params.label
% case 'spike'              % spike-threshold
%   switch params.technique
%     case 'threshold'
%       params.blockDur = 1;  % hours; amount of data to pull at once
%   end
% case 'burst'
%   switch params.technique
%     case 'linelength'     % burst-linelength
%       params.function = @(x) sum(abs(diff(x))); % sum(x.*x); % feature function
%       params.windowLength = 1;         % sec, duration of sliding window
%       params.windowDisplacement = 0.5;    % sec, amount to slide window
%       params.blockDurHr = 1;            % hours; amount of data to pull at once
%       params.smoothDur = 0;   % sec; width of smoothing window
%       params.minThresh = 2e5;    % X * stdev(signal); minimum threshold to detect burst; 
%       params.minDur = 2;    % sec; min duration of the seizures
%       params.addAnnotations = 1; % upload annotations to portal
%       params.viewData = 1;  % look at the data while it's running?
%       params.saveToDisk = 0;  % save calculations to disk?
%   end
% case 'seizure'
%   switch params.technique
%     case 'energy'     % seizure-area
%       params.function = @(x) sum(x.*x); % feature function
%       params.windowLength = 2;         % sec, duration of sliding window
%       params.windowDisplacement = 1;    % sec, amount to slide window
%       params.blockDurHr = 1;            % hours; amount of data to pull at once
%       params.smoothDur = 30;   % sec; width of smoothing window
%       params.minThresh = 15e8;    % X * stdev(signal); minimum threshold to detect burst; 
%       params.minDur = 10;    % sec; min duration of the seizures
%       params.addAnnotations = 1; % upload annotations to portal
%       params.viewData = 0;  % look at the data while it's running?
%       params.saveToDisk = 0;  % save calculations to disk?
%   end
% end


