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
labview.name     = 'Read LabView Files';
labview.val      = {filename,headerlength,timeind,decimalTriggerInds,binaryTriggerInds,trigThresh};
labview.help     = {'Reading LabView data'}';
labview.prog     = @lbv_read;
labview.vout     = @vout_lbv_read;
labview.modality = {'EEG'};
end
%==========================================================================
function labview = lbv_read(job)
% construct the S struct
S=job;
[a,b,~]=fileparts(S.filename{1});
outfile= fullfile(a,[ b, '.mat']);

data    = spm_opm_read_lvm(S);
save(outfile,'data');

labview.data= {outfile};


end
%==========================================================================
function dep = vout_lbv_read(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'Prepared Labview Data';
dep.src_output = substruct('.','data');
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end