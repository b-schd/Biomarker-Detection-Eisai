
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
    %{'EEG EEG1.3A-B','EEG EEG2.3A-B','EMG EMG.3'}, %not used
    {'EEG EEG1A-B','EEG EEG2A-B'},
    };

%anonymous functions
%EnergyFn = @(x) mean(x.^2);
%ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
LLFn = @(x,fs) nanmean(abs(diff(x)));

feature='freq';
switch feature
    case 'freq'
        featFn = @calc_featureswithfreqcorr;
        prefix = 'freq-global';
    case 'LL'
        featFn = LLFn;
        prefix = 'LL-global';
end
        

%% split layer based on channels
%[~,splitTimes,splitCh] = split_annotations(allEvents, timesUSec, eventChannels);
winLen = 2;
winDisp = 1;
mode = 'global'; %same # of electrodes, interchangeable
%% Train Model
% for each dataset
f_X = [];
f_Y = [];
for i = 1:numel(session.data)
    fprintf('Working on %s\n',session.data(i).snapName);
    feat = [];
    feat2 = [];
    fs = session.data(i).sampleRate;
    channels = session.data(i).channelLabels(:,1);
    layer_names = {session.data(i).annLayer.name};
    layer = layer_names(ismember(layer_names,'True_Seizures'));
    %if layer exists
    if ~isempty(layer)
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
    else
        fprintf('No Annotations\n');
    end
end
% model = TreeBagger(100,X,Y);
c = [0 50; 1 0];
model = fitcsvm(f_X,f_Y,'KernelFunction','linear','Cost',c);
%lr = mnrfit(X,categorical(Y+1))
cv = crossval(model);
kfoldLoss(cv)
%% detect for current dataset
for i = 1:numel(session.data)
    channels = session.data(i).channelLabels(:,1);
    %run through each mouse combination
    if mode=='global'
        for j = 1:numel(groupChannels)
            curCh = groupChannels{j};
            ch = find(ismember(channels,curCh));
            run_detections(session.data(i),model,winLen,winDisp,ch,featFn,prefix,'append')
        end
    end
end

%%
% Input: timesUsec, eventChannels
% output: cell array of timesUsec, eventChannels that correspond to
% unique channels
function [splitEvents, splitTimes, splitCh] = split_annotations(events,times,channels)   
    C =channels;
    maxLengthCell=max(cellfun('size',C,2));  %finding the longest vector in the cell array
    for i=1:length(C)
        for j=cellfun('size',C(i),2)+1:maxLengthCell
             C{i}(j)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
    A=cell2mat(C); %A is your matrix
    [~,~,IC] = unique(A,'rows','sorted');
    splitEvents = cell(1,max(IC));
    splitTimes = cell(1,max(IC));
    splitCh = cell(1,max(IC));
    for i=1:max(IC)
        splitEvents{i} = events(IC==i);
        splitTimes{i} = times(IC==i,:);
        tmp = cell2mat(channels(IC==i));
        splitCh{i} = tmp(1,:);
    end
end

function [feat, chs] = extractFeaturesFromAnnotationLayer(dataset,layerName,winLen,winDisp,fs,featFn)
    [~, timesUSec, chs] = getAnnotations(dataset,layerName);
    annotIdx = timesUSec/1e6*fs;
    %If duration of annotation is 0, assume 10 second duration
    %seizure
    duration= annotIdx(:,2)-annotIdx(:,1);
    annotIdx(duration==0,2) = annotIdx(duration==0,1)+10*fs;
    feat = cell(size(annotIdx,1),1);
    for k = 1:size(annotIdx,1)
        data = dataset.getvalues(annotIdx(k,1):annotIdx(k,2),chs{k});
        %extract features over windows
        %features = extract_features(data,fs)
        feat{k} = runFuncOnWin(data,fs,winLen,winDisp,featFn);
        %out{j} = cell2mat(runFuncOnWin(data,fs,2,1,@calc_featureswithfreqcorr));
    end
    % should all be the same since annotations should be split before
    % running
    ch = chs{1};
end

function run_detections(dataset,model,winLen,winDisp,ch,featFn,prefix,layerOption)
    datasetName = dataset.snapName;
    fs = dataset.sampleRate;

    startTime = 1; % in seconds
    stopTime = dataset.rawChannels(1).get_tsdetails.getDuration/1e6;
    
    params.timeOfInterest=[0 stopTime];
    params.filtFlag = 0; %1 to filter, 0 to ignore
    params.winLen = winLen; %(s)
    params.winDisp = winDisp; %(s)
    params.blockLen = 60*60*1; %Length of data to get from cloud at one time.
    params.IEEGid = 'hoameng';
    params.IEEGpwd = 'hoa_ieeglogin.bin';
%    outPrefix=calcFeature_v19_par(dataset,ch{1},params,features);
    %% extract features in segs
    blockLenSecs = 3600;
    numBlocks = ceil((stopTime-startTime)/blockLenSecs);
    startIdx = 1;
    szIdx = [];
    for i = 23:numBlocks
        startBlockPt = round(startIdx+(blockLenSecs*(i-1)*fs));
        endBlockPt = round(startIdx+blockLenSecs*i*fs-1);
        fprintf('Block %d of %d, ch%d_%d idx %d:%d\n',i,numBlocks,ch(1),ch(2),startBlockPt,endBlockPt);
        fsave = sprintf('%s-ch%d_%d-idx-%d-%d-%s.mat',dataset.snapName,ch(1),ch(2),startBlockPt,endBlockPt,prefix);
        if ~isempty(dir(fsave))
            fprintf('%s - mat found, loading\n',fsave);
            f = load(fsave);
            feat = f.feat;
        else
            data = dataset.getvalues(startBlockPt:endBlockPt,ch);
            feat = runFuncOnWin(data,fs,winLen,winDisp,featFn);
            save(fsave,'feat','-v7.3')
        end
        yhat = predict(model,feat);
        %yhat = cell2mat(yhat);
        idx = find(yhat==1);
        winIdx = idx*(winLen-(winLen-winDisp)); %now in secs
        tmpIdx = winIdx*fs + startBlockPt - 1;
        szIdx = [szIdx; tmpIdx];
    end
%     %% detect
%     fn = dir(outPrefix);
%     parFeats = load(fn.name);
%     parFeats = parFeats.parFeats;
%     %% store true detections
%     Yhat = predict(model,parFeats);
%     y = cell2mat(Yhat);
    
    if numel(szIdx)>0
        channels = cell(size(szIdx,1),1);
        for c = 1:numel(channels)
            channels{c} = ch;
        end
        szIdxRange = [szIdx szIdx+winLen*fs];
        %uploadAnnotations(dataset,sprintf('%s_detected_clips',prefix),szIdxRange/fs*1e6,channels,'SZ',layerOption)

        %duration features
        szIdx = sort(szIdx);
        dsz = diff(szIdx);
        finalSzIdx = [];
        tmpIdx = szIdx;
        sidx = 1;
        eidx = 1;
        for i = 1:numel(dsz)
            if dsz(i) > (winLen)*fs
                toAdd = [tmpIdx(sidx),tmpIdx(eidx)];
                finalSzIdx = [finalSzIdx; toAdd];
                sidx = eidx + 1;
                eidx = sidx;
            else
                eidx = eidx + 1;
            end
        end
        toAdd = [tmpIdx(sidx),tmpIdx(eidx)];
        finalSzIdx = [finalSzIdx; toAdd];

        durations = finalSzIdx(:,2)-finalSzIdx(:,1);
        finalSzIdx = finalSzIdx(durations>(winLen*fs*2),:);
        finalSzIdx(:,2) = finalSzIdx(:,2) + winLen; %correct for left shift detections

        channels = cell(size(finalSzIdx,1),1);
        for c = 1:numel(channels)
            channels{c} = ch;
        end
        uploadAnnotations(dataset,sprintf('%s_detected_seizures',prefix),finalSzIdx/fs*1e6,channels,'SZ',layerOption)
    end
end
