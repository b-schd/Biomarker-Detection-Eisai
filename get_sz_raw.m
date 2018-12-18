

addpath(genpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2'));
%addpath('~/gdriveshort/Libraries/Utilities/hline_vline');
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Analysis'))
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox/portalGit/Utilities'))
%javaaddpath('Z:\public\USERS\hoameng/Libraries/ieeg-matlab-1.13.2/IEEGToolbox/lib/ieeg-matlab.jar');

%params = initialize_task_humanNV;
params = initialize_task;


channelIdxs = cell(numel(session.data),1);
for i = 1:numel(session.data)
    channelIdxs{i} = [1 4];
end

blockLenSecs = 3600;
timeOfInterest = [];

% Load data
session = loadData(params);

dat = load('trainset.mat');
dat = dat.trainSet;
dat = table2array(dat);
dat =dat(dat(:,3)==1,:);
subj = unique(dat(:,1));
session = IEEGSession(sprintf('I033_A00%.2d_D001',subj(1)),'hoameng','hoa_ieeglogin.bin');
alldat = cell(numel(subj),1);
for i = 1:numel(subj)
    if i ~= 1
        session.openDataSet(sprintf('I033_A00%.2d_D001',subj(i)))
    end
    cursubj = subj(i);
    tmp = dat(dat(:,1)==cursubj(1,1),:); % find seizures
    if ~isempty(tmp)
        [~,timesUSec,eventChannels] = getAnnotations(session.data(i),'initialLL200-linelength');
        timesSec = timesUSec/1e6;
        subjdat = cell(size(tmp,2),1);
        fs =session.data(i).sampleRate;
        for j = 1:size(tmp,1)
            startSec = tmp(j,2);
            tmp1 = abs(timesSec(:,1)-startSec);
            idx = find(tmp1==min(tmp1));
            subjdat{j} = session.data(i).getvalues(timesUSec(idx,1)/1e6*fs:timesUSec(idx,2)/1e6*fs,[1 3]);
        end 
        alldat{i} = subjdat;
    end
end
    


dat = load('trainset2.mat');
dat = dat.tab;
dat = table2array(dat);
subj2 = unique(dat(:,1));
session = IEEGSession(sprintf('I033_A00%.2d_D001',subj2(1)),'hoameng','hoa_ieeglogin.bin');
alldat2 = cell(numel(subj2),1);
for i = 1:numel(subj2)
    if i~=1
        session.openDataSet(sprintf('I033_A00%.2d_D001',subj2(i)))
    end
    cursubj = subj2(i);
    tmp = dat(dat(:,1)==cursubj(1,1),:); % find seizures
    if ~isempty(tmp)
        subjdat = cell(size(tmp,2),1);
        fs =session.data(i).sampleRate;
        for j = 1:size(tmp,1)
            subjdat{j} = session.data(i).getvalues(tmp(j,2)*fs:tmp(j,3)*fs,[1 3]);
        end 
        alldat2{i} = subjdat;
    end
end

alldat = [alldat; alldat2];
subj = [subj; subj2];
save('allsz.mat','alldat','subj');


