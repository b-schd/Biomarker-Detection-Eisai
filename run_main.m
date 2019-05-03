%run_main
% Runs large seizure detector, 5 seconds in duration, 

prefix = 'test';
run_pipeline('kaggle','RFKaggle',prefix,25,5,1,1,'append')
