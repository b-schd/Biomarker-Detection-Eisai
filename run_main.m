%run_main
% Runs large seizure detector, 5 seconds in duration, 

prefix = 'test';

feature= 'kaggle';      % 'LL', 'freq', 'kaggle', or 'ts'
model = 'RFKaggle';     % 'SVM', 'RF', or 'RFKaggle'
layer_prefix = 'test'; 
inputC = 25;
durationThreshold = 5;  % seconds or 'min0.5' to set automatically 
indivMode = 1;
globalMode = 1;
layerOption = 'append';  % 'append' or 'overwrite'

run_pipeline(feature,'RFKaggle',layer_prefix,inputC,durationThreshold,indivMode,globalMode,layerOption)

% run generate_table.m