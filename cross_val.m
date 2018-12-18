load('ii_featswithcorr.mat');
load('sz_featswithcorr.mat');

%% CROSS VALIDATION

subj = unique(subjIdx);
sz = unique(szIdx);
indices = crossvalind('Kfold',numel(unique(szIdx)),numel(unique(szIdx)));
ii_idx = 1:72%numel(sz);
ii_idx = repmat(ii_idx,60*2/2,1);
ii_idx = reshape(ii_idx,numel(ii_idx),1);
%ii_idx 
ii_idx_group = crossvalind('Kfold',1:72,numel(sz));
for i = 1:max(unique(indices))
    testSZ = szIdx(indices==i); 
    Xtest_seiz = sz_feats(szIdx==testSZ,:);
    Xtest_ii = ii_feats(ii_idx==ii_idx_group(i),:);
    Xtest = [Xtest_seiz;Xtest_ii];
    Ytest = [ones(size(Xtest_seiz,1),1);zeros(size(Xtest_ii,1),1)];
    
    trainSZ = szIdx(indices~=i); 
    Xtrain_seiz = sz_feats(szIdx~=testSZ,:);
    Xtrain_ii = ii_feats(ii_idx~=ii_idx_group(i),:);
    Xtrain = [Xtrain_seiz;Xtrain_ii];
    Ytrain = [ones(size(Xtrain_seiz,1),1);zeros(size(Xtrain_ii,1),1)];
    
    %train model
    rfmodel = TreeBagger(300,Xtrain,Ytrain,...
    'Method','classification','oobvarimp','on','Cost',[0 1; 1 0]);
    
    %test model
    [Ypred Yscore] = predict(rfmodel,Xtest);
    
    %post_processing
    C = confusionmat(categorical(Ytest),categorical(Ypred))

    
    
end



figure
plot(oobError(rfmodel));
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Classification Error');

vars = cell(18,1);
k = 1;
for ch = 1:2
    for i = 1:18
        vars{k} = sprintf('%d-%d',ch,i);
        k = k+1;
    end
end
varimp = rfmodel.OOBPermutedVarDeltaError';
[~,idxvarimp]= sort(varimp);
labels = vars(idxvarimp);

figure
barh(varimp(idxvarimp),1); ylim([1 52]);
set(gca, 'YTickLabel',labels, 'YTick',1:numel(labels))
title('Variable Importance'); xlabel('score')





