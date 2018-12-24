
function run_detections(dataset,model,winLen,winDisp,durations,ch,featFn,prefix,layerOption)
    datasetName = dataset.snapName;
    fs = dataset.sampleRate;

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
    for i = 1:numBlocks
        startBlockPt = round(startIdx+(blockLenSecs*(i-1)*fs));
        endBlockPt = round(startIdx+blockLenSecs*i*fs-1);
        fprintf('Block %d of %d, ch%d_%d idx %d:%d...',i,numBlocks,ch(1),ch(2),startBlockPt,endBlockPt);
        fsave = sprintf('%s-%s,ch%d_%d-idx-%d-%d.mat',dataset.snapName,prefix,ch(1),ch(2),startBlockPt,endBlockPt);
        if ~isempty(dir(fsave))
            fprintf('%s - mat found, loading\n',fsave);
            f = load(fsave);
            feat = f.feat;
        else
            fprintf('...running\n')
            data = dataset.getvalues(startBlockPt:endBlockPt,ch);
            feat = runFuncOnWin(data,fs,winLen,winDisp,featFn);
            save(fsave,'feat','-v7.3')
        end
        yhat = predict(model,feat);
        %yhat = cell2mat(yhat);
        idx = find(yhat==1);
        winIdx = idx*(winLen-(winLen-winDisp)); %now in secs
        tmpIdx = winIdx*fs + startBlockPt - 1;
        szIdx = [szIdx; tmpIdx];
    end
%     %% detect
%     fn = dir(outPrefix);
%     parFeats = load(fn.name);
%     parFeats = parFeats.parFeats;
%     %% store true detections
%     Yhat = predict(model,parFeats);
%     y = cell2mat(Yhat);
    
    if numel(szIdx)>0
        channels = cell(size(szIdx,1),1);
        for c = 1:numel(channels)
            channels{c} = ch;
        end
        szIdxRange = [szIdx szIdx+winLen*fs];
        %uploadAnnotations(dataset,sprintf('%s_detected_clips',prefix),szIdxRange/fs*1e6,channels,'SZ',layerOption)

        %duration features
        szIdx = sort(szIdx);
        dsz = diff(szIdx);
        finalSzIdx = [];
        tmpIdx = szIdx;
        sidx = 1;
        eidx = 1;
        for i = 1:numel(dsz)
            if dsz(i) > (winLen)*fs
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
        finalSzIdx = finalSzIdx(detdurations>(winLen*fs*2),:); %duration threshold >2xwinLen
        detdurations = finalSzIdx(:,2)-finalSzIdx(:,1);
        finalSzIdx = finalSzIdx(detdurations>(mean(durations)*0.95*fs),:); %duration threshold min of training durations
        fprintf('Duration threshold: %d seconds\n',mean(durations)*0.95)
        if size(finalSzIdx,1) > 0

            channels = cell(size(finalSzIdx,1),1);
            for c = 1:numel(channels)
                channels{c} = ch;
            end
            uploadAnnotations(dataset,sprintf('%s_detected_seizures',prefix),finalSzIdx/fs*1e6,channels,'SZ',layerOption)
        else
            fprintf('Detections removed based on duration threshold\n')
        end
    end
end