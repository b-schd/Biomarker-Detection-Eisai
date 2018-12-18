dat = load('allsz.mat');

subj = dat.subj;
dat = dat.alldat;
allout = cell(numel(subj),1);
winLen = 2;
winDisp = 2;
fs = 2000;
        
subjIdx = [];
szIdx = [];
% for each subj
k = 1;
for i = 1:numel(subj)
    %for each sz
    tmp = dat{i};
    out = cell(numel(tmp),1);
    for j = 1:numel(tmp)
        tmpdat = tmp{j};
        
        %calculate features
        %out = [out cell2mat(runFuncOnWin(dat,fs,winLen,winDisp,@featFn))];
        out{j} = cell2mat(runFuncOnWin(tmpdat,fs,winLen,winDisp,@calc_featureswithfreqcorr));
        subjIdx = [subjIdx ;ones(size(out{j},1),1)*subj(i)];
        szIdx = [szIdx;ones(size(out{j},1),1)*k];
        k = k +1;
    end
    allout{i} = cell2mat(out);
end

sz_feats = cell2mat(allout);
save('sz_featswithfreqcorr.mat','sz_feats','subjIdx','szIdx');

dat = load('interictal.mat');
dat = dat.interictals;
out = cell(numel(dat),1);
for j = 1:numel(dat)
    fprintf('%d of %d\n',j,numel(dat))
    tmpdat = dat{j};

    winLen = 2;
    winDisp = 2;
    fs = 2000;
    %calculate features
    %out = [out cell2mat(runFuncOnWin(dat,fs,winLen,winDisp,@featFn))];
    out{j} = cell2mat(runFuncOnWin(tmpdat,fs,winLen,winDisp,@calc_featureswithfreqcorr));

end
ii_feats = cellfun(@(x)cell2mat(x),out,'UniformOutput',0);
ii_feats = cell2mat(ii_feats);
save('ii_featswithfreqcorr.mat','ii_feats');


