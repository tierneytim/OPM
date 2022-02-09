function [pD] = spm_opm_bandpow(S)
% computed a sliding window of power in  particular frquency band 
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                 - Default: no Default
%   S.win           - size of window in ms            - Default: no Default              
%   S.band          - band to compute power in        - Default: no Default              
% Output:
%   pD               - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id$

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),      error('an MEEG object must be supplied.'); end
if ~isfield(S, 'win'),    error('a time window(ms) must be supplied.'); end
if ~isfield(S, 'band'),   error('a frequency band must be supplied'); end


%-Filter to specific band
%--------------------------------------------------------------------------
args = [];
args.D = S.D;
args.type = 'butterworth';
args.band = 'bandpass';
args.freq = S.band;
args.dir = 'twopass';
args.order = 2;

fD = spm_eeg_filter(args);


%-get number of samples 
%--------------------------------------------------------------------------
nsamples = S.D.fsample*S.win/1000;
if mod(nsamples,2)==0
    nsamples = nsamples +1;
end



%-compute power 
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),['pow_' fname(S.D)]);
pD = clone(S.D,outname);
pD.save();

ind = indchantype(fD,'MEG');
otherind = setdiff(1:size(fD,1),ind);

pow = movvar(squeeze(fD(ind,:,1))',nsamples);

pD(ind,:,:)= sqrt(pow');
pD(otherind,:,:)= S.D(otherind,:,:);
pD.save();

end 
