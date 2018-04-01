%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)
%% Sensor Level Example - data supplied 
S =[];
S.data = ones(10,2)*5;
S.fs =1000;
D = spm_opm_create(S);

%% Sensor Level Example - data supplied with labels
S =[];
S.data = ones(80,2)*5;
S.fs =1000;
S.pinout = 'OPMpinout_20171018';
S.sensorsUsed = 'OPM2cast_MedianNerve';
D = spm_opm_create(S);

%% Sensor Level Example - data supplied with labels and triggers
S =[];
S.data = ones(80,2)*5;
S.trig = ones(3,2)*1;
S.fs =1000;
S.pinout = 'OPMpinout_20171018';
S.sensorsUsed = 'OPM2cast_MedianNerve';
D = spm_opm_create(S);

%% Source Level Example - labels, triggers, positions and MRI file
S =[];
S.data = ones(80,2)*5;
S.trig = ones(3,2)*1;
S.fs =1000;
S.pinout = 'OPMpinout_20171018';
S.sensorsUsed = 'OPM2cast_MedianNerve';
S.pos = 'SEF_coarse';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% Simulation Example - no data, but positions and MRI file supplied
S =[];
S.pos = 'SEF_coarse';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% Simulation Example - Individual MRI, custom spacing between sensors
S =[];
S.space = 15;
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% Simulation Example - Individual MRI, custom mesh file
S =[];
S.space = 15;
S.cortex='testCustom.gii';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% Simulation Example - Individual MRI, wholehead coverage
S =[];
S.space = 15;
S.wholehead = 1;
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% Simulation Example - Template MRI file, custom  sensor spacing
S =[];
S.space = 15;
D = spm_opm_create(S);

%% Simulation Example - Template MRI file, custom mesh
S =[];
S.cortex = 'testCustom.gii';
D = spm_opm_create(S);

%% Simulation Example - Template MRI file, custom offsets (SQUID vs OPM) 
%results consistent with Iivanainen et al (2017) and  Boto et al (2017) 

S =[];
S.space = 30;
S.offset = 1; 
S.lead = 1;
[D,L1] = spm_opm_create(S);
sc1 = D.inv{1}.forward.scale;
OPM1 = mean(max(abs(L1)))/sc1;

S =[];
S.space = 30;
S.offset = 20; 
S.lead = 1;
[D,L2] = spm_opm_create(S);
sc2 = D.inv{1}.forward.scale;
squid1 = mean(max(abs(L2)))/sc2;


S =[];
S.space = 30;
S.offset = 6; 
S.lead = 1;
[D,L1] = spm_opm_create(S);
sc1 = D.inv{1}.forward.scale;
OPM2 = mean(max(abs(L1)))/sc1;

S =[];
S.space = 30;
S.offset = 50; 
S.lead = 1;
[D,L2] = spm_opm_create(S);
sc2 = D.inv{1}.forward.scale;
squid2 = mean(max(abs(L2)))/sc2;