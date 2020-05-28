# Biomarker Detection

Code to automatically run seizure detection based on training annotation layers. This was developed with a project on seizure detection in a collaboration with the Jensen lab at UPenn and Eisai, and later in collaboration with the Wolf lab.

Currently this project is aimed at seizure detection, but will be expanded for other biomarker detection tasks in the future.

## Overview of Algorithm
The algorithm extracts features from labeled seizures and non-seizures and uses several algorithms to make detections on the entire recording. Seizures and non-seizures were manually annotated by animal. Each labeled annotation is broken down into 2 second window lengths with 1 second displacements, and feature sets are extracted from all channels. Thus, each example consisted of features extracted from two seconds of EEG from allchannels. These feature sets varied from classic time domain features such as line length to power spectra and correlations in the time and frequency domain. Features were then extracted from the entire recordings each animal at a time and predictions were made based on a chosen classifier. Classifiers included random forests and support vector machines, with parameters manually tuned for specificity. Since each classification is two seconds in duration, a duration threshold was set to capture major seizures of interest. At the time of the initial classification, mean - 2 SD of true seizure durations was roughly 12 seconds, and the initial algorithm was thus ran with a duration threshold of 12 seconds. Predictions were made based on an individual animal (Each prediction only made from annotations from the same animal) as well as globally (model was trained on all animals). The cross-validation accuracy for the best algorithm (based on the Kaggle competition) was roughly 91% individually by window, 100% PPV by seizure.

## Pitfalls
1. Due to differences in activity/recording quality between rats, individually trained algorithms performed better than global algorithms, as expected.  
2. The kaggle algorithm leverages correlations. Because this is a two channel setup per animal, normalizations of the power spectra were done differently than originally documented (within channel instead of across channel). 
3. Based on training/test accuracies, more data is required to improve the global classifier (versus more complex classifier). At a certain point, manual review may be more efficient.
4. The initial algorithm's marked durations are not accurate, as the algorithm was tuned to detect the number of seizures, not necessarily duration. Proceed with caution when interpreting durations.

## Requirements: 
ieeg-matlab toolbox (https://bitbucket.org/ieeg/ieeg/downloads/)  
portal-matlab-tools (https://github.com/ieeg-portal/portal-matlab-tools.git)  
Biomarker-Detection-IEEG (https://github.com/b-schd/Biomarker-Detection-IEEG)  

## Setup:
Ensure that all packages/repos above are in one directory along with extracted matlab toolbox.  
(e.g. Projects/Biomarker-Detection-IEEG, Projects/portal-matlab-tools, Projects/ieeg-matlab-1.xx.x)  
1. git clone https://github.com/b-schd/Biomarker-Detection-IEEG 
2. git clone https://github.com/ieeg-portal/portal-matlab-tools.git  
3. Download ieeg-matlab toolbxo from link above  
4. Setup ieeg-matlab toolbox. (see portal-matlab-tools/IEEGTutorial.m). Once you set up your account, the login.bin needs to be in the matlab path. It may be easiest to place it in the ieeg-matlab folder.

## Run locally on matlab  
The function below uses relative paths to the other packages. Navigate to the Eisai folder within matlab and run the function below:    
1. Edit initialize_task.m to include relevant parameters, datasets, grouped channels 
2. run_pipeline('kaggle','RFKaggle','test',25,12,1,1,'append') (run kaggle algorithm, prefix layernames with 'test', set C parameter to 25, minimum duration of 12 seconds, run individually as well as global, and append to existing annotations if any)  


# Misc tips
1. removeAnnotations from portal-matlab-tools may be useful in removing layers that have been added.  
2. To calculate duration of seizures,  
a. session = IEEGSession('[dataset]','username','password')  
b. annotations = getAnnotations(IEEGDataset object, 'layer name')  
c. ([annotations.end] - [annotations.start])/1e6 %durations of all annotations in seconds (note that these durations are expected to be shorter than the duration of the actual seizures)  



