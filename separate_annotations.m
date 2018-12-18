
addpath(genpath('../ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('../portal-matlab-tools/Analysis'))
addpath(genpath('../portal-matlab-tools/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

params = initialize_task;

% Load data
session = loadData(params);

% 
for j = 1:numel(session.data)
    [allEvents, timesUSec, channels] = getAnnotations(session.data(j),'EDF Annotations');
    startTimes = [];
    endTimes = [];
    eventChannels = [];

    for i = 1:numel(allEvents)
        if strcmp(allEvents(i).description,'Seizure')
            startTimes = [startTimes timesUSec(i,1)];
            endTimes = [endTimes timesUSec(i,2)];
            eventChannels = [eventChannels channels(i)];
        end
    end
    if ~isempty(startTimes)
        uploadAnnotations(session.data(j),'True_Seizures', [startTimes' endTimes'], eventChannels, 'Seizure','overwrite')
    end
end