fn = dir('Jensen_Eisai*.mat');

for i =1:numel(fn)
    curfn = fn(i).name;
    t = load(curfn);
    feat = t.feat;
    curfn = strrep(curfn,'groupA-ch1-5-idx','groupA-ch1_5-idx');
    fprintf('%s\n',curfn)
    save(curfn,'feat')
end