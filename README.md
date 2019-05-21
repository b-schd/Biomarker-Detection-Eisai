# Eisai

Code to automatically run seizure detection based on training annotation layers

## Requirements: 
ieeg-matlab toolbox (https://bitbucket.org/ieeg/ieeg/downloads/)  
portal-matlab-tools (https://github.com/ieeg-portal/portal-matlab-tools.git)  
Eisai.git (https://github.com/hoameng/Eisai.git)  

## Setup:
Ensure that all packages/repos above are in one directory along with extracted matlab toolbox.  
(e.g. Projects/Eisai, Projects/portal-matlab-tools, Projects/ieeg-matlab-1.xx.x)  
1. git clone https://github.com/hoameng/Eisai.git  
2. git clone https://github.com/ieeg-portal/portal-matlab-tools.git  
3. Download ieeg-matlab toolbxo from link above  
4. Setup ieeg-matlab toolbox. (see portal-matlab-tools/IEEGTutorial.m). Once you set up your account, the login.bin needs to be in the matlab path. It may be easiest to place it in the ieeg-matlab folder.

## Run locally on matlab  
1. run_pipeline('kaggle','RFKaggle','test',25,12,1,1,'append') (run kaggle algorithm, prefix layers of test, set C parameter to 25, minimum duration of 12 seconds, run individually as well as global, and append to existing annotations if any)

## Run on linux:
Note, these scripts were tailored for a specific project in mind. To do: generalize
1. Edit initialize_task.m to include relevant parameters, datasets, grouped channels, as well as 
2. run_main - run a "large seizure" detector  
a. open ssh client (putty/cygwin/terminal)  
b. type: ssh username@borel.seas.upenn.edu, enter your password  
c. navigate to your Eisai folder   
d. matlab -nodisplay    
e. run_pipeline('kaggle','RFKaggle','test',25,12,1,1,'append')  

tips: google "linux screen" 

# Misc tips
1. removeAnnotations from portal-matlab-tools may be useful in removing layers that have been added.  
2. To calculate duration of seizures,  
a. session = IEEGSession('[dataset]','username','password')  
b. annotations = getAnnotations(IEEGDataset object, 'layer name')  
c. ([annotations.end] - [annotations.start])/1e6 %durations of all annotations in seconds (note that these durations are expected to be shorter than the duration of the actual seizures)  


