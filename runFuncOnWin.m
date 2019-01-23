function out = runFuncOnWin(x,fs,winLen,winDisp,featFn)
% Usage:  out = runFuncOnWin(x,fs,winLen,winDisp,featFn)
% Function function featFn on moving windows.
% Input:
%   x       :   N x Ch
%   fs      :   sample Rate
%   winLen  :   window length in seconds
%   winDisp :   window displacement in seconds

%anonymous functions
%EnergyFn = @(x) mean(x.^2);
%ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
%LLFn = @(x,fs) nanmean(abs(diff(x)));
NumWins = @(xLen, fs, winLen, winDisp)floor((xLen-(winLen-winDisp)*fs)/(winDisp*fs));

numWindows = NumWins(length(x),fs,winLen,winDisp);

out = cell(numWindows,1);
maxLength = 0;
nanFlag = 0;
for i = 1:numWindows
    %find current window
    data = x(1+(i-1)*(winDisp*fs):((i-1)*(winDisp)+winLen)*fs,:);
    try
        out{i} = featFn(data,fs);
        if maxLength == 0
            maxLength = numel(out{i});
        end
    catch e
        fprintf('Unable to calculate for window %d\n',i);
        fprintf(1,'Identifier:\n%s',e.identifier);
        fprintf(1,'Error! The message was:\n%s',e.message);
        out{i} = NaN;
        nanFlag = 1;
    end
end
if nanFlag
    for i = 1:numel(out)
        if isnan(out{i})
            out{i} = NaN(1,maxLength);
        end
    end
end
out = cell2mat(out);
