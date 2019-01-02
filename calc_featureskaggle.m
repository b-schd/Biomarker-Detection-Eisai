
function feats = calc_featureskaggle(data, fs)
%needs generalization to more than 2 channels

%     %power from 1-48hz 
%     freqRange = 1:48;
%     P = pmtm(data,[],size(data,1),fs);
%     feats = log10(abs(P(freqRange,:)))';
%     feats = reshape(feats, 1,[]);
   
    P = pmtm(data,[],size(data,1),fs);
    P2 = P(1:47,:);
    
    feats = log10(abs(P2));
    feats = reshape(feats, 1,[]);
    
    %normalize across frequency buckets per channel (*note, different
    %than what is described in MHills documentation. Necessary due to
    %2 channels, where he described normalizing across channels would
    %result in correlation of 1 due to only 2 channels
    P3 = (P2 - mean(P2,1))./std(P2,0,1);
    r = corr(P3);
    e = abs(eig(r));
    r = r(1,2);
    
    %correlation between channels
    data2 = (data - mean(data,1))./std(data,0,1);
    
    r2 = corr(data2);
    e2 = abs(eig(r2));
    r2 = r2(1,2);
    
    feats = [feats e' r e2' r];
    
end