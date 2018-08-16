%% Housekeeping
clear all
addpath('spm12')
spm('defaults', 'eeg')
addpath('OPM\')
dir = 'OPM\testData';
cd(dir)
%% read lvm
S=[];
S.filename='QZFM_test.zip';
S.nchannels=81;
S.timeind=1;
S.decimalTriggerInds=74:80;
S.binaryTriggerInds=81;
S.trigThresh=4;
lbv=spm_opm_read_lvm(S);

%% read Opti
S=[];
S.fname='optiTest.csv';
S.optiTrigger=lbv.binaryTrigs;
S.triggerTime=lbv.time;
opti = spm_opm_read_optitrack(S);

%%
S =[];
S.data = lbv.B';
S.trig = lbv.decimalTrigs';
S.fs =1200;
S.scale = 1e6/2.7;
S.pinout= 'OPMpinout_20171018';
S.sensorsUsed='OPM_2_cast_opti_test';
S.other=opti.motion';
S.pos = 'wholeheadcast_opti_test';
S.sMRI= 'msMQ0484_orig.img';
D = spm_opm_create(S);

%% filter the data
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [1 80];
S.dir = 'twopass';
S.order = 2;
D = spm_eeg_filter(S);

%% epoch the data
S =[];
S.D=D;
S.timewin=[-100 300];
D= spm_opm_epochTrigger(S);

%% denoising
S=[];
S.D=D;
S.confounds={'REF';'PHYS'};
S.gs=0;
S.derivative=1;
dD = spm_opm_synth_gradiometer(S);

