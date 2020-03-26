% upload_test_layer

% Read in test_annots.csv
annots_file='../test_annots.csv';

datasetName= params.datasetID{1};
session= IEEGSession(datasetName, 'bscheid', '/Users/bscheid/ieeg_pwd.bin');

anns= readtable(annots_file);
anns_ieeg = IEEGAnnotation.createAnnotations([anns.start]*1e6, [anns.stop]*1e6, anns.event)
newLayer = session.data.addAnnLayer(params.marked_seizure_layer);
newLayer.add(anns_ieeg);