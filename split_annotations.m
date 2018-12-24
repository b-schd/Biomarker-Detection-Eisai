
%%
% Input: timesUsec, eventChannels
% output: cell array of timesUsec, eventChannels that correspond to
% unique channels
function [splitEvents, splitTimes, splitCh] = split_annotations(events,times,channels)   
    C =channels;
    maxLengthCell=max(cellfun('size',C,2));  %finding the longest vector in the cell array
    for i=1:length(C)
        for j=cellfun('size',C(i),2)+1:maxLengthCell
             C{i}(j)=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
    A=cell2mat(C); %A is your matrix
    [~,~,IC] = unique(A,'rows','sorted');
    splitEvents = cell(1,max(IC));
    splitTimes = cell(1,max(IC));
    splitCh = cell(1,max(IC));
    for i=1:max(IC)
        splitEvents{i} = events(IC==i);
        splitTimes{i} = times(IC==i,:);
        tmp = cell2mat(channels(IC==i));
        splitCh{i} = tmp(1,:);
    end
end