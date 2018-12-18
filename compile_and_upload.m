
addpath(genpath('../../../Libraries/ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('portal-matlab-tools/Analysis'))
addpath(genpath('portal-matlab-tools/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

%params = initialize_task_humanNV;
params = initialize_task;

% Load data
session = loadData(params);

channelIdxs = cell(numel(session.data),1);
for i = 1:numel(session.data)
    channelIdxs{i} = [1 3];
end

startTime = 1;
numSZ = [];
for subjIdx = 1:numel(session.data)
    dataset = session.data(subjIdx);
    datasetName = session.data(subjIdx).snapName;
    %load all parFeats
    fn = dir(sprintf('tmp2/%s-scores-par*',datasetName));
    fn = {fn.name};
    preds = [];
    for fnIdx = 1:numel(fn)
        tmpparfeat = load(sprintf('tmp2/%s',fn{fnIdx}));
        preds= [preds; tmpparfeat.parFeats];
    end
   % final_preds = cell2mat(preds);
    temp = conv(preds,ones(1,5),'same');
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
            timer = 180;
        end
    end
    
    szTimes = szTimes'*1e6;
    eventChannels = {};
    for i = 1:numel(szTimes)
        eventChannels{i} = [1 3];
    end
    uploadAnnotations(session.data(subjIdx),'Detected SZ-FC',szTimes,eventChannels,'SZ',1)
    toAdd = [subjIdx size(szTimes,1)];
    numSZ = [numSZ; toAdd];
end

