
function feats = calc_features(data, fs)

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
    
    
    %correlation between channels
    %feats = [feats triu(corrcoef(data),-1)];
    
    %cross correlation
    [acor, lag] = xcorr(data,.25*fs,'coeff');
    lag = lag(acor(:,2)==max(acor(:,2)));
    first_lag = lag(1);
    feats = [feats max(acor(:,2)) first_lag];
    
end