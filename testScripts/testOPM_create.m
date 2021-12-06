%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)
%% Sensor Level Example - data supplied 
S =[];
S.data = 'meg.bin';
D = spm_opm_create(S);
%% Sensor Example - should be same as previous
S =[];
S.data = 'meg.bin';
S.coordystem='coordsystem.json';
S.positions='positions.tsv';
S.channels='channels.tsv';
S.mef='meg.json';
D = spm_opm_create(S);
%% merge example
S =[];
S.data = 'meg.bin';
S.sMRI= 'T1w.nii';
D1 = spm_opm_create(S);

S =[];
S.data = 'meg.bin';
S.sMRI= 'T1w.nii';
D2 = spm_opm_create(S);


fis = [fullfile(D1);fullfile(D2)];
S= [];
S.D = fis;
S.recode= 'same';
D3 = spm_eeg_merge(S);
%% Source Example 
S =[];
S.data = 'meg.bin';
S.sMRI='T1w.nii';
D = spm_opm_create(S);

%% Source example with no structural MRI, but 3 fiducials specified in coordsystem.json
S               = [];
S.data          = 'meg.bin';
S.coordystem    = 'coordsystem.json';
S.positions     = 'positions.tsv';
S.channels      = 'channels.tsv';
S.template      = 1;
S.meshres       = 3;
D               = spm_opm_create(S);

%% Simulation Example - default
S =[];
D = spm_opm_create(S);

%% Simulation Example - no data, but positions and MRI file supplied
S =[];
S.positions = 'positions.tsv';
S.sMRI= 'T1w.nii';
D = spm_opm_sim(S);

%% Simulation Example - Individual MRI, custom spacing between sensors
S =[];
S.space = 25;
S.sMRI= 'T1w.nii';
D = spm_opm_sim(S);

%% Simulation Example - Individual MRI, custom mesh file
S =[];
S.space = 25;
S.cortex='testCustom.gii';
S.sMRI= 'T1w.nii';
D = spm_opm_sim(S);

%% Simulation Example - Individual MRI, not wholehead
S =[];
S.space = 25;
S.wholehead = 0;
S.sMRI= 'T1w.nii';
D = spm_opm_sim(S);

%% Simulation Example - Template MRI file, custom  sensor spacing
S =[];
S.space = 25;
D = spm_opm_sim(S);

%% Simulation Example - Template MRI file, custom mesh
S =[];
S.cortex = 'testCustom.gii';
D = spm_opm_sim(S);

%% Simulation Example - individual subject,fixed positions
S =[];
S.positions = 'positions.tsv';
S.sMRI= 'T1w.nii';
D = spm_opm_sim(S);

%% Fieldtrip to SPM Example
load('data_ft.mat');

S = [];
S.data = data_ft;
S.sMRI='T1w.nii';
D = spm_opm_create_ft(S);

