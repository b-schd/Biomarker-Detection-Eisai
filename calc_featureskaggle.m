
function feats = calc_featureskaggle(data,varargin)
%data should be NxC, where rows are samples and columns are channels

%     %power from 1-48hz 
%     freqRange = 1:48;
%     P = pmtm(data,[],size(data,1),fs);
%     feats = log10(abs(P(freqRange,:)))';
%     feats = reshape(feats, 1,[]);

    data1 = TimeCorrelation(data);
    data2 = FreqCorrelation(data);
   

    feats = [data1;data2]';
    
end

function out = TimeCorrelation(data)
    %correlation between channels, normalize if > 2 channels
    if size(data,2) > 2
        data = (data - mean(data,2))./std(data,1,2);
    end
    %center to mean and scale by unit variance
    r2 = corr(data);
    e2 = abs(eig(r2));
    ref = r2*0+1;
    t = triu(r2,1);
    ref = triu(ref,1);
    t = t(ref==1);
    out = [e2;t];
end

function out = FreqCorrelation(data)
    P = rfft(data);
    P2 = P(2:48,:); %skip 0 hz
    P2 = abs(P2);
    P2 = log10(P2);
    
    %normalize across frequency buckets per channel (*note, different
    %than what is described in MHills documentation. Necessary due to
    %2 channels, where he described normalizing across channels would
    %result in correlation of 1 due to only 2 channels. I believe his
    %documentation is incorrect
    P3 = P2;
    if size(P3,2)>2
        P3 = (P3 - mean(P3,2))./std(P3,1,2);
    end
    r2 = corr(P3);
    e2 = abs(eig(r2));
    ref = r2*0+1;
    t = triu(r2,1);
    ref = triu(ref,1);
    t = t(ref==1);
    
    P2 = reshape(P2,[],1);
    out = [P2;e2;t];
end
function rfft = rfft(a)
     ffta = fft(a);
     rfft = ffta(1:(floor(length(ffta)/2)+1),:);
end