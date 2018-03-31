function labview = spm_cfg_opm_read_lvm
% configuration file for reading lab view file
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Tim Tierney

%--------------------------------------------------------------------------
% labview file
%--------------------------------------------------------------------------
filename        = cfg_files;
filename.tag    = 'filename';
filename.name   = 'File Name';
filename.filter = '(.lvm|.zip)';
filename.num    = [1 1];
filename.help   = {'Select the (zipped) lvm file.'};

%--------------------------------------------------------------------------
% nchannels
%--------------------------------------------------------------------------
nchannels       = cfg_entry;
nchannels.tag     = 'nchannels';
nchannels.name    = 'No. of channels';
nchannels.help    = {'The number of columns containing data in the labview file.'};
nchannels.strtype = 'r';
nchannels.num     = [1,1];
nchannels.val     = {81};

%--------------------------------------------------------------------------
% headerlength
%--------------------------------------------------------------------------
headerlength         = cfg_entry;
headerlength.tag     = 'headerlength';
headerlength.name    = 'No. of header lines';
headerlength.help    = {'The number of lines of text containing header information'};
headerlength.strtype = 'r';
headerlength.num     = [1,1];
headerlength.val     = {23};



%--------------------------------------------------------------------------
% timeind
%--------------------------------------------------------------------------
timeind         = cfg_entry;
timeind.tag     = 'timeind';
timeind.name    = 'Time Index';
timeind.help    = {'The column number of the timevariable'};
timeind.strtype = 'r';
timeind.num     = [1,1];
timeind.val     = {1};

%--------------------------------------------------------------------------
% Decimal Triggers
%--------------------------------------------------------------------------
decimalTriggerInds         = cfg_entry;
decimalTriggerInds.tag     = 'dec';
decimalTriggerInds.name    = 'Index of Decimal Triggers';
decimalTriggerInds.help    = {'Column numbers of channels that should be interpreted as decimal triggers'};
decimalTriggerInds.strtype = 'r';
decimalTriggerInds.num     = [1,Inf];
decimalTriggerInds.val     = {74:81};

%--------------------------------------------------------------------------
% Binary Triggers
%--------------------------------------------------------------------------
binaryTriggerInds         = cfg_entry;
binaryTriggerInds.tag     = 'bin';
binaryTriggerInds.name    = 'Index of Decimal Triggers';
binaryTriggerInds.help    = {'Column numbers of channels that should be interpreted as decimal triggers'};
binaryTriggerInds.strtype = 'r';
binaryTriggerInds.num     = [1,Inf];
binaryTriggerInds.val     = {74:81};

%--------------------------------------------------------------------------
% Trigger Threshld
%--------------------------------------------------------------------------
trigThresh         = cfg_entry;
trigThresh.tag     = 'thresh';
trigThresh.name    = 'Trigger Threshold';
trigThresh.help    = {'Threshold to apply to triggers.'};
trigThresh.strtype = 'r';
trigThresh.num     = [1,1];
trigThresh.val     = {4};

%--------------------------------------------------------------------------
% read
%--------------------------------------------------------------------------
labview          = cfg_exbranch;
labview.tag      = 'labview';
labview.name     = 'LabView';
labview.val      = {filename, nchannels,headerlength,timeind,decimalTriggerInds,binaryTriggerInds,trigThresh};
labview.help     = {'Reading LabView data'}';
labview.prog     = @lbv_read;
labview.vout     = @vout_lbv_read;
labview.modality = {'EEG'};

%==========================================================================
function labview = lbv_read(job)
% construct the S struct
S=job;
labview.lbv    = spm_opm_read_lvm(S);

%==========================================================================
function dep = vout_lbv_read(job)
% return dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'LabView object';
dep(1).src_output = substruct('.','lbv');
dep(1).tgt_spec   = cfg_findspec({{'strtype','lbv'}});
