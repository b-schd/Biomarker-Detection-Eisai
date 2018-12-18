fn = dir('Jensen_Eisai*0.mat');

for i =1:numel(fn)
    curfn = fn(i).name;
    t = load(curfn);
    feat = t.feat;
    curfn = strrep(curfn,'0.mat','0-LL.mat');
    fprintf('%s\n',curfn)
    save(curfn,'feat')
end