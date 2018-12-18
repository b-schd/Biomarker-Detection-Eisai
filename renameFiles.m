dataset = [2 5 7 11 12 15 35 36 37 38 66 67 68 69 70 76]
for dataset = 76
    fn = dir(sprintf('tmp2/I033_A00%0.2d_D001-scores-par*',dataset));
    fn = {fn.name};

    for i = 1:numel(fn)
        tmp = load(sprintf('tmp2/%s',fn{i}));
        id = str2num(fn{i}(27:end-4));
        if id < 1000
            parFeats = tmp.parFeats;
            newfn = sprintf('tmp2/I033_A00%0.2d_D001-scores-par%0.4d.mat',dataset,id);
            save(newfn,'parFeats');
        end
        delete(sprintf('tmp2/%s',fn{i}));
    end
end