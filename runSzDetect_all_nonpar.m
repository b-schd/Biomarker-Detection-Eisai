
addpath(genpath('../../../Libraries/ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('portal-matlab-tools/Analysis'))
addpath(genpath('portal-matlab-tools/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

%params = initialize_task_humanNV;
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
try
 %parpool(20);
catch
end
load 'szMdl.mat';
for subjIdx = 1:numel(session.data)
    dataset = session.data(subjIdx);
    datasetName = session.data(subjIdx).snapName;
    fs = dataset.sampleRate;

    startTime = 1; % in seconds
    stopTime = session.data(subjIdx).rawChannels(1).get_tsdetails.getDuration/1e6;

    numSegs = ceil((stopTime-startTime)/(2000));
    start = startTime;
    preds = cell(numSegs,1); % will be a list of classification scores for each 2-second window
    session2 = IEEGSession(datasetName,'hoameng','hoa_ieeglogin.bin');
    dataset = session2.data;
    for i = 1:numSegs
        fprintf('Working on segment %g of %g...\n',i,numSegs);
       
        % pull the next 2000 seconds
        start = startTime + (i-1)*2000;
        endTime = min(start+2000,stopTime);
        pulled=false;
        while pulled==false
            try
                curData = dataset.getvalues(start*fs:endTime*fs-1,[1 3]);
                pulled = true;
            catch
                fprintf('Pulled failed... trying again...\n');
            end
        end

        % break into 2-second clips

        numClips = floor(size(curData,1)/(2*fs));
        tempStart = 1;
        tmpfeats = {};
        for j = 1:numClips
        % skip clips with all zero lead or NaNs (but keep track so the times
        % line up
            clip = curData(tempStart:tempStart+2*fs-1,:);
            if any(any(isnan(clip)))
                tmpfeats{j} = NaN(1,18);
                continue 
            end
            if any(all(clip==0))
                tmpfeats{j} = NaN(1,18);
                continue
            end
            tmpfeats{j} = calc_features(clip,fs);
            tempStart = tempStart+2*fs;
        end
        feats = cell2mat(tmpfeats');
        [~, tmpscore] = predict(mdl,feats);
        [r c] = find(isnan(feats));
        tmpscore(r,:) = NaN;
        preds{i} = tmpscore(:,2);
    end
    save(sprintf('%s-predsnonpar',datasetName),'preds');
    final_preds = cell2mat(preds);
    temp = conv(final_preds,ones(1,5),'same');
    finalPreds = temp>2;
    %%
    szTimes = [];
    timer = 0;
    for i = 1:length(finalPreds)
        if timer > 0
            timer = timer - 2;
            continue
        end
        if finalPreds(i)
            szTimes = [szTimes startTime+2*i];
            timer = 300;
        end
    end
    
    szTimes = szTimes'*1e6;
    eventChannels = {};
    for i = 1:numel(szTimes)
        eventChannels{i} = [1 3];
    end
    uploadAnnotations(session.data(subjIdx),'Detected SZ',szTimes,eventChannels,'SZ',0)
end

