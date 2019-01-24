
function run(feature,model,layer_prefix, durationThreshold, indivMode,globalMode,layerOption)
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
        model= 'RF';
        feature_prefix = 'freq';
    case 'kaggle'
        featFn = @calc_featureskaggle;
        model = 'RFkaggle';
        feature_prefix = 'kgl';
    case 'LL'
        featFn = LLFn;
        model = 'SVM';
        feature_prefix = 'LL';
    case 'ts'
        featFn = @calc_tfeatures;
        model = 'SVM';
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
for i = 1:numel(session.data)
    fprintf('Working on %s\n',session.data(i).snapName);
    feat = [];
    feat2 = [];
    tmpszdurations = [];
    fs = session.data(i).sampleRate;
    channels = session.data(i).channelLabels(:,1);
    layer_names = {session.data(i).annLayer.name};
    layer = layer_names(ismember(layer_names,'True_Seizures'));
    %if layer exists
    if ~isempty(layer)
        [~, timesUSec, chs] = getAnnotations(session.data(i),layer);
        tmp = timesUSec(:,2)-timesUSec(:,1);
        tmpszdurations = [tmpszdurations; tmp/1e6];
        [feat, ch] = extractFeaturesFromAnnotationLayer(session.data(i),layer{1},winLen,winDisp,fs,featFn);
    end
    layer = layer_names(ismember(layer_names,'Non_Seizures'));
    %if layer exists
    if ~isempty(layer)
        [feat2, ch2] = extractFeaturesFromAnnotationLayer(session.data(i),layer{1},winLen,winDisp,fs,featFn);
    end
    
    if ~isempty(feat) && ~isempty(feat2)
        %train model
        feat = cell2mat(feat);
        feat2 =cell2mat(feat2);
        X = [feat; feat2];
        Y = [ones(size(feat,1),1); zeros(size(feat2,1),1)];
        f_X = [f_X;X];
        f_Y = [f_Y;Y];
        
        if indivMode==1
            switch model
                case 'SVM'
                    C = [0 50; 1 0];
                    mdl = fitcsvm(X,Y,'KernelFunction','linear','Cost',C);
                case 'RF'
                    mdl = TreeBagger(1000,X,Y,'OOBPrediction','off');
                case 'RFkaggle'
                    C = [0 numel(Y)/(numel(unique(Y))*sum(Y==1)); numel(Y)/(numel(unique(Y))*sum(Y==0)) 0];
                    % n_samples / (n_classes * np.bincount(y))
                    mdl = TreeBagger(3000,X,Y,'MinLeafSize',2,'Cost',C,'SampleWithReplacement','on','OOBPrediction','on');
                    %kaggle solution below
                    %mdl = TreeBagger(3000,X,Y,'MinLeafSize',2,'SampleWithReplacement','off','OOBPrediction','off');
            end
            %lr = mnrfit(X,categorical(Y+1))
            %cv = crossval(model);
            %kfoldLoss(cv)
            %% detect for current dataset
            if strcmp(durationThreshold,'min0.5') || isempty(durationThreshold)
                durationThreshold = min(tmpszdurations)*0.5;
            end
            run_detections(session.data(i),mdl,winLen,winDisp,durationThreshold,minThreshold,mergeThreshold,ch{1},featFn,layer_prefix,feature_prefix,'indiv',layerOption)
        end
    else
        fprintf('No Annotations\n');
    end
    szdurations = [szdurations; tmpszdurations];
end
switch model
    case 'SVM'
        c = [0 50; 1 0];
        mdl = fitcsvm(f_X,f_Y,'KernelFunction','linear','Cost',c);
    case 'RF'
        mdl = Treebagger(1000,f_X,f_Y,'OOBPrediction','off');
    case 'RFkaggle'
        C = [0 numel(Y)/(numel(unique(Y))*sum(Y==1)); numel(Y)/(numel(unique(Y))*sum(Y==0)) 0];
        % n_samples / (n_classes * np.bincount(y))
        mdl = TreeBagger(3000,X,Y,'MinLeafSize',2,'Cost',C,'SampleWithReplacement','on','OOBPrediction','on');
        %mdl = TreeBagger(3000,X,Y,'MinLeafSize',2,'SampleWithReplacement','off','OOBPrediction','off');
end
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
                run_detections(session.data(i),mdl,winLen,winDisp,durationThreshold, minThreshold,mergeThreshold,ch,featFn,layer_prefix,feature_prefix,'global',layerOption)
            end
        end
    end
end



end
