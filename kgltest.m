data = csvread('testszdata.csv');
truefeats = csvread('testszkglfeat.csv');
fs = 400


feats = calc_featureskaggle(data, fs)