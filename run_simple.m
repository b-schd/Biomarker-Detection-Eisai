
function run_simple(feature,layer_prefix, durationThreshold, indivMode,globalMode,layerOption)
%
% Usage: run(feature,model,layer_prefix, durationThreshold, indivMode,globalMode,layerOption)
% durationThreshold : int in seconds or 'min0.5' to set automatically
% to half of annotated detections., default is 'min0.5'

%feature == 'LL' or 'freq'
%model = 'SVM' or 'RF'
%layerOption = 'append' or 'overwrite'
%default: run('LL','SVM',1,1,'append')
addpath(genpath('../ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('../portal-matlab-tools/Analysis'))
addpath(genpath('../portal-matlab-tools/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

params = initialize_task;

% Load data
session = loadData(params);
% 
% GENERATE TABLE
%find length of each dataset
subj = cell(numel(session.data),1);
dR = zeros(numel(session.data),1);
for i = 1:numel(session.data)
    subj{i} = session.data(i).snapName;
    dR(i) = session.data(i).rawChannels(1).get_tsdetails.getDuration/1e6/60/60/24;
end

channelIdxs = cell(numel(session.data),1);
for i = 1:numel(session.data)
    channelIdxs{i} = [1 3];
end

% groupChannels = {
%     {'EEG EEG1.1A-B','EEG EEG2.1A-B','EMG EMG.1'},  
%     {'EEG EEG1.2A-B','EEG EEG2.2A-B','EMG EMG.2'},
%     %{'EEG EEG1.3A-B','EEG EEG2.3A-B','EMG EMG.3'}, %not used
%     {'EEG EEG1A-B','EEG EEG2A-B','EMG EMG'},
%     };

groupChannels = {
    {'EEG EEG1.1A-B','EEG EEG2.1A-B'},  
    {'EEG EEG1.2A-B','EEG EEG2.2A-B'},
    {'EEG EEG1.3A-B','EEG EEG2.3A-B'}, %not used
    {'EEG EEG1A-B','EEG EEG2A-B'},
    };

%anonymous functions
%EnergyFn = @(x) mean(x.^2);
%ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
LLFn = @(x,fs) nanmean(abs(diff(x)));

% PARAMETERS
%feature='LL';
%indivMode = 1;
%globalMode = 1;
winLen = 2;
winDisp = 1;
mergeThreshold = winLen; %in seconds was 2*winLen
minThreshold = winLen;
switch feature
    case 'freq'
        featFn = @calc_featureswithfreqcorr;
       % model= 'RF';
        feature_prefix = 'freq';
    case 'kaggle'
        featFn = @calc_featureskaggle;
        %model = 'RFkaggle';
        feature_prefix = 'kgl';
    case 'LL'
        featFn = LLFn;
        %model = 'SVM';
        feature_prefix = 'LL';
    case 'ts'
        featFn = @calc_tfeatures;
        %model = 'SVM';
        feature_prefix = 'ts';
end
        

%% split layer based on channels
%for i = 1:numel(session.data)
%    splitAnnotationsByChannel(session.data(i),'True_Seizures');
%end

%% Train Model
% for each dataset
f_X = [];
f_Y = [];
szdurations = [];

%% detect for current dataset
for i = 1:numel(session.data)
    channels = session.data(i).channelLabels(:,1);
    %run through each mouse combination
    if globalMode==1
        for j = 1:numel(groupChannels)
            curCh = groupChannels{j}; %get channels for each mouse
            ch = find(ismember(channels,curCh))'; %find channel numbers
            if ~isempty(ch)
                if strcmp(durationThreshold,'min0.5') || isempty(durationThreshold)
                    durationThreshold = min(szdurations)*0.5;
                end
                run_detections_simple(session.data(i),winLen,winDisp,durationThreshold, minThreshold,mergeThreshold,ch,featFn,layer_prefix,feature_prefix,'global',layerOption)
            end
        end
    end
end



end
