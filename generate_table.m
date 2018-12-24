%% generate table 


% Group ChannelGroup SeizureStart SeizureEnd Duration

%feature == 'LL' or 'freq'
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
collate = [];
for i = 1:numel(session.data)
    subj = session.data(i).snapName;
    layer_names = {session.data(i).annLayer.name};
    idx = ~cellfun(@isempty,regexpi(layer_names, 'seizures$'));
    layers = layer_names(idx);
    for j = 1:numel(layers)
        [~, timesUsec, ch] = getAnnotations(session.data(i),layers{j});
        timesUsec = round(timesUsec/1e6,3);
        for k = 1:size(timesUsec,1)
            duration = (timesUsec(k,2)- timesUsec(k,1));
            toAdd = {subj ch{k} layers{j} timesUsec(k,1) timesUsec(k,2) duration};
            collate = [collate; toAdd];
        end
    end
end
T = cell2table(collate);
T.Properties.VariableNames = {'Group' 'Ch_num1' 'Ch_num2' 'Layer' 'Start' 'End' 'Duration'};
writetable(T,'table.csv')