function splitAnnotationsByChannel(dataset,layer_name)
%% 
% Function will find layer_name in dataset and output layer_name_x, where x = 1, 2, ...
% for each combination of channels associated with each annotation
%%

[annots, timesUSec, channels] = getAnnotations(dataset,layer_name);

existingCh = {};
newTimesUSec = {};
newCh = {};
for i = 1:numel(channels)
    ch = num2str(channels{i});
    idx = find(ismember(ch,existingCh));
    if ~isempty(idx)
        newTimesUSec{idx} = [newTimesUSec{idx}; timesUSec(i,:)];
        newCh{idx} = [newCh{idx} channels(i)];
    else
        existingCh = [existingCh ch];
        newTimesUSec = [newTimesUSec timesUSec(i,:)];
        newCh = [newCh channels(i)];
    end
end

if size(newTimesUSec,1) > 1
    %more than 1 different annotation
    fprintf('%d unique channels found in %s layer in %s..splitting',size(newTimesUSec,1), layer_name,dataset.snapname)
    for i = 1:numel(newTimesUSec)
        uploadAnnotations(dataset,sprintf('%s_%d',layer_name,i),newTimesUSec{i},newCh{i},annots(1).type,'overwrite')
    end
else
    fprintf('No distinct channels found in %s\n',dataset.snapname)
end
