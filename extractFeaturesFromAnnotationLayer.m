function [feat, chs] = extractFeaturesFromAnnotationLayer(dataset,layerName,winLen,winDisp,fs,featFn)
    [~, timesUSec, chs] = getAnnotations(dataset,layerName);
    annotIdx = timesUSec/1e6*fs;
    %If duration of annotation is 0, assume 10 second duration
    %seizure
    duration= annotIdx(:,2)-annotIdx(:,1);
    annotIdx(duration==0,2) = annotIdx(duration==0,1)+10*fs;
    feat = cell(size(annotIdx,1),1);
    for k = 1:size(annotIdx,1)
        data = dataset.getvalues(annotIdx(k,1):annotIdx(k,2),chs{k});
        %extract features over windows
        %features = extract_features(data,fs)
        feat{k} = runFuncOnWin(data,fs,winLen,winDisp,featFn);
        %out{j} = cell2mat(runFuncOnWin(data,fs,2,1,@calc_featureswithfreqcorr));
    end
    % should all be the same since annotations should be split before
    % running
    ch = chs{1};
end