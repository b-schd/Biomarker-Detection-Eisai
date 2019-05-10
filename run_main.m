%run_main

minDuration = 12;
C = 25;

prefix = sprintf('%dsec-c%d,minDuration,C);
indivMode = 1;
globalMode = 1;

run_pipeline('kaggle','RFKaggle',prefix,C,minDuration,indivMode,globalMode,'append')

% kaggle is the algorithm used in the kaggle competition 
% RFKaggle are the features and model used in this algorithm (random forest)
% prefix is the prefix of the layer
% C = parameter that controls false positive rate (higher = lower false positives, lower sensitivity)
% minDuration = minimum duration of a seizure (default 12)
% indivMode = 1
% globalMode = 1 
