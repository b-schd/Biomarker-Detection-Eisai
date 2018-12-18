function out = run_calc_features_all(subnum)
%% Establish IEEG Sessions
% Establish IEEG Sessions through the IEEGPortal. This will allow on demand
% data access
addpath(genpath('~/gdriveshort/Libraries/ieeg-matlab-1.13.2-borel'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('~/gdriveshort/Projects/p05-IEEGPortalToolbox/portalGit/Analysis'))
addpath(genpath('~/gdriveshort/Projects/p05-IEEGPortalToolbox/portalGit/Utilities'))
%javaaddpath('../../Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

addpath(genpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Analysis'))
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

%params = initialize_task_humanNV;
params = initialize_task_NV_2;

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
% % 
% % %create seizure mat file
% layername{1} = 'Seizure_CCS';
% layername{2} = 'Seizures';
% times = cell(numel(session),1);
% annTimes = [];
% for i = 1:numel(session.data)
%     for ln = 1:numel(layername)
%         try
%            [~, times] = getAnnotations(session.data(i),layername{ln}); 
%            times = ceil(times/1e6/60/60/24);
%            toAdd = [ones(size(times,1),1)*i times];
%            annTimes = [annTimes;toAdd];
%         catch
%         end
%     end
% end



channelIdxs = cell(numel(session.data),1);
for i = 1:numel(session.data)
    channelIdxs{i} = 1:16;
end

params.timeOfInterest=[];%[0 60*60*24];
params.filt.order = 3;
params.filt.wn = 190;
params.filt.type ='low';
params.filtFlag = 1;
params.winLen = 60*5;
params.winDisp = 60*5;
params.blockLen = 15*60*1; 
for i = subnum
    %calcFeature_v14_par(session.data(i),channelIdxs{i},params,{'LL','hwamp','area','power','energy','amp'});s
    %calcFeature_v18_par(session.data(i),channelIdxs{i},params,{'LL','rms','hw','area','power','energy'});
    calcFeature_v18_par(session.data(i),channelIdxs{i},params,{'hw'});
end
%job = batch('batchscript',0,{session.data,params,'power'},'Pool',1);
