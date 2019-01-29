
function run_detections_simple(dataset,winLen,winDisp,durationThreshold,minThreshold,mergeThreshold,ch,featFn,layer_prefix,prefix,prefix2,layerOption)
%prefix 1 for feature name, prefix 2 for global/indiv (to reuse .mat
%feature calculations and to differentiate layers
    bg_win = 180;
    std_thres = 5;
    datasetName = dataset.snapName;
    fs = dataset.sampleRate;
    base_path = 'tmp/';
    startTime = 1; % in seconds
    stopTime = dataset.rawChannels(1).get_tsdetails.getDuration/1e6;
    
    params.timeOfInterest=[0 stopTime];
    params.filtFlag = 0; %1 to filter, 0 to ignore
    params.winLen = winLen; %(s)
    params.winDisp = winDisp; %(s)
    params.blockLen = 60*60*1; %Length of data to get from cloud at one time.
    params.IEEGid = 'hoameng';
    params.IEEGpwd = 'hoa_ieeglogin.bin';
%    outPrefix=calcFeature_v19_par(dataset,ch{1},params,features);
    %% extract features in segs
    blockLenSecs = 3600;
    numBlocks = ceil((stopTime-startTime)/blockLenSecs);
    startIdx = 1;
    szIdx = [];
    yhat = [];
    cur_features = [];
    for i = 1:numBlocks
        startBlockPt = round(startIdx+(blockLenSecs*(i-1)*fs));
        endBlockPt = round(startIdx+blockLenSecs*i*fs-1);
        fprintf('Block %d of %d, ch%d_%d idx %d:%d...',i,numBlocks,ch(1),ch(2),startBlockPt,endBlockPt);
        fprintf('...running\n')
        data = dataset.getvalues(startBlockPt:endBlockPt,ch);
        feat = runFuncOnWin(data,fs,winLen,winDisp,featFn);
        for win = 1:size(feat,1)
            yhat = [yhat; 0];
            if any(isnan(data))
                continue
            else
                if size(cur_features,1)<bg_win
                    cur_features = [cur_features; feat(win)];
                else
                    meanfeat = mean(cur_features);
                    sdfeat = std(cur_features);
                    if (feat(win) > meanfeat+std_thres*sdfeat)
                        yhat(end) = 1;
                        fprintf('Possible seizure\n');
                    end
                    cur_features = [cur_features(2:end,:);feat(win)];
                end
            end
        end
    end
    
    idx = find(yhat==1);
    winIdx = idx*(winLen-(winLen-winDisp)); %now in secs
    tmpIdx = winIdx*fs + startBlockPt - 1;
    szIdx = [szIdx; tmpIdx];
  
%     %% detect
%     fn = dir(outPrefix);
%     parFeats = load(fn.name);
%     parFeats = parFeats.parFeats;
%     %% store true detections
%     Yhat = predict(model,parFeats);
%     y = cell2mat(Yhat);
    
    if numel(szIdx)>0
%         channels = cell(size(szIdx,1),1);
%         for c = 1:numel(channels)
%             channels{c} = ch;
%         end
        szIdxRange = [szIdx szIdx+winLen*fs]; %set duration
        %uploadAnnotations(dataset,sprintf('%s_detected_clips',prefix),szIdxRange/fs*1e6,channels,'SZ',layerOption)

        %combine detections less than 1 second apart
        szIdx = sort(szIdx);
        dsz = diff(szIdx);
        finalSzIdx = [];
        tmpIdx = szIdx;
        sidx = 1;
        eidx = 1;
        for i = 1:numel(dsz)
            if dsz(i) > mergeThreshold*fs
                toAdd = [tmpIdx(sidx),tmpIdx(eidx)];
                finalSzIdx = [finalSzIdx; toAdd];
                sidx = eidx + 1;
                eidx = sidx;
            else
                eidx = eidx + 1;
            end
        end
        toAdd = [tmpIdx(sidx),tmpIdx(eidx)];
        finalSzIdx = [finalSzIdx; toAdd];

        finalSzIdx(:,2) = finalSzIdx(:,2) + winLen*fs; %correct for left shift detections

        detdurations = finalSzIdx(:,2)-finalSzIdx(:,1);
        finalSzIdx = finalSzIdx(detdurations>minThreshold * fs,:); %min threshold
        detdurations = finalSzIdx(:,2)-finalSzIdx(:,1);
        finalSzIdx = finalSzIdx(detdurations>(durationThreshold*fs),:); %duration threshold min of training durations
        fprintf('Duration threshold: %d seconds\n',durationThreshold)
        if size(finalSzIdx,1) > 0
            channels = cell(size(finalSzIdx,1),1);
            for c = 1:numel(channels)
                channels{c} = ch;
            end
            uploadAnnotations(dataset,sprintf('%s-%s-%s_detected',layer_prefix,prefix,prefix2),finalSzIdx/fs*1e6,channels,'SZ',layerOption)
        else
            fprintf('Detections removed based on thresholds\n')
        end
    end
end