load 'szMdl.mat';
subjNum = 36;
sessionName = ['I033_A00' num2str(subjNum) '_D001'];
session = IEEGSession(sessionName, 'sbaldassano', 'sba_ieeglogin.bin');
dataset = session.data;
fs = dataset.sampleRate;

startTime = 190793; % in seconds
stopTime = 208390;

numSegs = ceil((stopTime-startTime)/(2000));
start = startTime;
preds = []; % will be a list of classification scores for each 2-second window
for i = 1:numSegs
    fprintf('Working on segment %g of %g...',i,numSegs);
    tic
    % pull the next 2000 seconds
    endTime = min(start+2000,stopTime);
    pulled=false;
    while pulled==false
        try
            curData = dataset.getvalues(start*fs:endTime*fs-1,[1 3]);
            pulled = true;
        catch
            fprintf('Pulled failed... trying again...\n');
        end
    end
    
    % break into 2-second clips
    
    numClips = floor(size(curData,1)/(2*fs));
    tempStart = 1;
    feats = [];
    for j = 1:numClips
    % skip clips with all zero lead or NaNs (but keep track so the times
    % line up
        clip = curData(tempStart:tempStart+2*fs-1,:);
        if any(any(isnan(clip)))
            preds = [preds; NaN];
            continue 
        end
        if any(all(clip==0))
            preds = [preds; NaN];
            continue
        end
        feats = [feats; calc_features(clip,fs)];
        tempStart = tempStart+2*fs;
    end
    [pred, score] = predict(mdl,feats);
    preds = [preds; score(:,2)];
    start=start+2000;
    toc
end
        
temp = conv(preds,ones(1,5),'same');
finalPreds = temp>2;
%%
szTimes = [];
timer = 0;
for i = 1:length(finalPreds)
    if timer > 0
        timer = timer - 2;
        continue
    end
    if finalPreds(i)
        szTimes = [szTimes startTime+2*i];
        timer = 300;
    end
end

