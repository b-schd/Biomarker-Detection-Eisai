
function feats = calc_tfeatures(data, fs)

    AreaFn = @(x) nanmean(abs(x));
    EnergyFn = @(x) nanmean(x.^2);
    ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
    LLFn = @(x,fs) nanmean(abs(diff(x)));

    feats = [AreaFn(data) EnergyFn(data) ZCFn(data) LLFn(data,fs)];
    
end