function create = spm_cfg_opm_create
% configuration file for creating OPM objects
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% Tim Tierney

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
data        = cfg_files;
data.tag    = 'data';
data.name   = 'File name for data';
data.filter = '.*.mat$';
data.num    = [1 1];
data.help   = {'Select the labview.mat file.'};

%--------------------------------------------------------------------------
% Sampling Frequency
%--------------------------------------------------------------------------
fs         = cfg_entry;
fs.tag     = 'fs';
fs.name    = 'Sampling Frequency';
fs.help    = {'The sampling frequency of the data in Hz'};
fs.strtype = 'r';
fs.num     = [1,1];
fs.val     = {1200};
%--------------------------------------------------------------------------
% Sampling Frequency
%--------------------------------------------------------------------------
scale        = cfg_entry;
scale.tag     = 'scale';
scale.name    = 'Scale Factor';
scale.help    = {'Scale factor (multiplied by data)to convert to units of fT'};
scale.strtype = 'r';
scale.num     = [1,1];
scale.val     = {1e6/2.7};

%--------------------------------------------------------------------------
% create
%--------------------------------------------------------------------------
create          = cfg_exbranch;
create.tag      = 'create';
create.name     = 'Create OPM object';
create.val      = {data,fs,scale};
create.help     = {'Create/simulate OPM data'}';
create.prog     = @opm_create;
create.vout     = @vout_opm_create;
create.modality = {'EEG'};
end
%==========================================================================
function out = opm_create(job)
% construct the S struct
S=job;

if(strmatch(job.data,''))
   S= rmfield(job,'data');
else
    load(job.data{1});
    S.data= data.B';
    trigs = (size(data.decimalTrigs,2)+size(data.binaryTrigs,2))>0;
    if(trigs)
    S.trigs= [data.decimalTrigs,data.binaryTrigs]';
    end
end

[a,b,c]=fileparts(job.data{1});
outfile= fullfile(a,['SPM_',b,'.dat']);
S.fname = outfile;

out.D= spm_opm_create(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};
end
%==========================================================================
function dep = vout_opm_create(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'OPM Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'OPM Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
end