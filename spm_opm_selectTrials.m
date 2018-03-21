function D = spm_opm_selectTrials(S)
%--------------------------------------------------------------------------
% Selects trials and places them in a subsetted MEG object
%--------------------------------------------------------------------------
% FORMAT D = spm_opm_selectTrials(S)
%   S               - input structure
% Optional fields of S:
%   S.D             - SPM meeg object    - Default: no default
%   S.inds          - trials to subset   - Default: bootstraps the trials
%   S.fname         - output file name   - Default: prefix with 'sub_'
% Output:
%  D           - MEEG object (also written to disk)
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%-Initialise nad argument check
%--------------------------------------------------------------------------

% check D supplied
if ~isfield(S, 'D')
    error('MEEG object needs to be supplied'); 
end

% if no inds supplied bootstrap
if ~isfield(S, 'inds')
    trials = 1:ntrials(S.D);
    S.inds   = datasample(trials,ntrials(S.D));
end

% if no filename supplied make one up
if ~isfield(S, 'fname')
    [a,b,c] = fileparts(fullfile(S.D));
    S.fname   = fullfile(a,['sub_',b,c]); 
end

%-get new dimensions and Condition Labels
%--------------------------------------------------------------------------
origDim = size(S.D);
dim = origDim;
dim(end)= length(S.inds);
condLabels = conditions(S.D,S.inds);

%-clone dataset and reset condition labels
%--------------------------------------------------------------------------
D = clone(S.D,S.fname,dim,2);
if(length(size(D))==3)
D(:,:,:)=S.D(:,:,S.inds);
elseif(length(size(D))==4)
D(:,:,:,:)=S.D(:,:,:,S.inds);
end
D = conditions(D,1:ntrials(D),condLabels);
save(D)
%--------------------------------------------------------------------------

end



