

addpath(genpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Analysis'))
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

%params = initialize_task_humanNV;
params = initialize_task;

session = loadData(params);
channelIdxs = cell(numel(session.data),1);

for i = 1:numel(session.data)
    channelIdxs{i} = [1 3];
end

%% Initialization
params.burst.blockLenSecs = 7200;
params.timeOfInterest = [];
params.IEEGid = 'hoameng';
params.IEEGpwd = 'hoa_ieeglogin.bin';
params.burst.minDur = .5;
params.burst.maxDur = 30;
params.burst.winLen = 1;
params.burst.thres =2;
params.burst.maxThres = 100;

% Load data
session = loadData(params);

for i = 1:numel(session.data)
     fprintf('Detecting in : %s\n',session.data(i).snapName);
     try
        filtFlag = 0;
        [burstTimes, burstChannels] = burst_detector_v3(session.data(i),channelIdxs{i},params);
        
        %uploadAnnotations(session.data(i), 'burst_detections',burstTimes,burstChannels,'burst');
     catch ME
        disp(ME)
        fprintf('Failed detecting in: %s\n',session.data(i).snapName);
     end
end

%% CLUSTER BURSTS ONCE DETECTED ABOVE
%INITIAL CLUSTER
idxByDataset = clusterAllBursts(session.data,'burst_detections');
for i = 1:numel(session.data)
    idx = idxByDataset{i}; %GET IDXS FOR DATASET
   	annots = getAllAnnots(session.data(i),'burst_detections'); %GET ALL ANNOTS
    for j=1:max(idx) %FOR EACH CLUSTER
        fprintf('Adding cluster %d...',j)
        newLayerName = sprintf('burst_detections_%d',j);
        try
            session.data(i).removeAnnLayer(newLayerName); %TRY TO REMOVE IN CASE ALREADY PRESENT
        catch
        end
        annLayer = session.data(i).addAnnLayer(newLayerName);
        ann=annots(idx==j);
        numAnnot = numel(ann);
        startIdx = 1;
        %add annotations 5000 at a time (freezes if adding too many)
        for k = 1:ceil(numAnnot/5000)
            fprintf('Adding %d to %d\n',startIdx,min(startIdx+5000,numAnnot));
            annLayer.add(ann(startIdx:min(startIdx+5000,numAnnot)));
            startIdx = startIdx+5000;
        end
        fprintf('...done!\n')
    end
end