
function feats = calc_tfeatures(data, fs)

AreaFn = @(x) nanmean(abs(x));
EnergyFn = @(x) nanmean(x.^2);
ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
LLFn = @
   
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