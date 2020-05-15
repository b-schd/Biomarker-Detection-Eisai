% upload_test_layer

% Read in annotations from spreadsheet
params = initialize_task;

datasetName= params.datasetID{1};
session= IEEGSession(datasetName, params.IEEGid, params.IEEGpwd);

anns= readtable(params.annotFile);

% Get annotations for class 1 (seizures)
i_anns=~isnan([anns.sz_start]); 
anns_ieeg{1} = IEEGAnnotation.createAnnotations([anns.sz_start(i_anns)]*1e6,...
    [anns.sz_stop(i_anns)]*1e6, anns.sz_event(i_anns));

% Get annotations for class 2 (usu. interictal periods)
i_anns=~isnan([anns.class2_start]);
anns_ieeg{2} = IEEGAnnotation.createAnnotations([anns.class2_start(i_anns)]*1e6,...
    [anns.class2_stop(i_anns)]*1e6, anns.class2_event(i_anns));

for i_layer=1:length(anns_ieeg)
    
    try
        newLayer = session.data.addAnnLayer(params.marked_seizure_layers{i_layer});
    catch ME
        
        if strcmp(ME.message, sprintf('Layer with name ''%s'' already exists.', params.marked_seizure_layers{i_layer}))
            ipt= input(sprintf('Layer %s already exists, append annotations to layer (a) or overwrite layer (o)?\n',...
              params.marked_seizure_layers{i_layer}), 's');
          
            if strcmp(ipt,'a')
                layernames= {session.data.annLayer.name};
                newLayer=session.data.annLayer(strcmp(layernames, params.marked_seizure_layers{i_layer})); 
                
            elseif strcmp(ipt, 'o')
                session.data.removeAnnLayer(params.marked_seizure_layers{i_layer});
                newLayer = session.data.addAnnLayer(params.marked_seizure_layers{i_layer}); 
                
            end
            
        else, throw(ME)
            
        end

    end
    newLayer.add(anns_ieeg{i_layer});
    fprintf('%d annotations uploaded to %s \n', length(anns_ieeg{i_layer}), params.marked_seizure_layers{i_layer})

end
