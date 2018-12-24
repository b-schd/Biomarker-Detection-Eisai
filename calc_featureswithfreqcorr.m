
function feats = calc_featureswithfreqcorr(data, fs)

%     %power from 1-48hz 
%     freqRange = 1:48;
%     P = pmtm(data,[],size(data,1),fs);
%     feats = log10(abs(P(freqRange,:)))';
%     feats = reshape(feats, 1,[]);
   
    P = pmtm(data,[],size(data,1),fs);
    
    binSize = 5; %hz
    freqBins = 1:binSize:size(P,1);
    Pbins = zeros(numel(freqBins),2);
    for i = 1:numel(freqBins)-1
        Pbins(i,:) = sum(P(freqBins(i):freqBins(i+1),:));
    end
    
    Pbins = log10(abs(Pbins(1:8,:)));
    feats = reshape(Pbins, 1,[]);
    r = corr(P);
    r = r(1,2);
    
    
    %correlation between channels
    %feats = [feats triu(corrcoef(data),-1)];
    
    %cross correlation
    [acor, lag] = xcorr(data,.25*fs,'coeff');
    lags = lag(acor(:,2)==max(acor(:,2)));
    feats = [feats max(acor(:,2)) lags(1) r];
    
end