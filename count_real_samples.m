
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

%need to clear 7
for subjIdx = 1:numel(session.data)
    dataset = session.data(subjIdx);
    datasetName = session.data(subjIdx).snapName;
    fs = dataset.sampleRate;

    startTime = 1; % in seconds
    stopTime = session.data(subjIdx).rawChannels(1).get_tsdetails.getDuration/1e6;
    if strcmp(dataset.snapName,'I033_A0005_D001')
        stopTime = 370000;
    end

    numSegs = ceil((stopTime-startTime)/(2000));
    start = startTime;
    numSamples = 0;
    for i = 1:numSegs
            fprintf('Working on segment %g of %g...\n',i,numSegs);

            % pull the next 2000 seconds
            start = startTime + (i-1)*2000;
            endTime = min(start+2000,stopTime);

            count = 0;
            pulled = false;
            curdata = 0;
            while count<10 && ~pulled
                try
                    curData = dataset.getvalues(start*fs:endTime*fs-1,[1 3]);
                    pulled = true;
                catch
                    fprintf('Pulled failed... trying again...\n');
                    count = count + 1;
                end
            end
            a = curData(:,1);
            numSamples = numSamples + sum(~isnan(a));
            
            % break into 2-second clips
    end
    snapName = dataset.snapName;
    save(sprintf('%s_numsample',snapName),'snapName','numSamples');
end

