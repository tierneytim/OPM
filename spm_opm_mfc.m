function [mfD,Yinds] = spm_opm_mfc(S)
% remove interference that behaves as if it was from a mean (magnetic) field
% FORMAT D = spm_opm_mfc(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.usebadchans   - logical to use channels marked as bad          - Default: 0
% Output:
%   D               - denoised MEEG object (also written to disk)
%   Yinds           - the indices of filtered channels
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_mfc.m 7646 2019-07-25 13:58:46Z tim $

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'usebadchans'),   S.usebadchans = 0; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'badChanThresh'), S.badChanThresh = 50; end
if ~isfield(S, 'balance'),       S.balance = 0; end

%-Get design matrix
%--------------------------------------------------------------------------
try
    s= sensors(S.D,'MEG');
catch
    error('Could not find sensor positions')
end

if(S.usebadchans)
    usedLabs= s.label;
else
    badLabels = chanlabels(S.D,badchannels(S.D));
    indsRem = [];
    for i =1:length(badLabels)
        indsRem=[indsRem strmatch(badLabels{i},s.label)];
    end
    LabInds = 1:length(s.label);
    sinds = setdiff(LabInds,indsRem);
    usedLabs= s.label(sinds);
end
X= s.coilori(sinds,:);
fprintf('%-40s: %30s\n','Created Design Matrix',spm('time'));

%-Compute projector
%--------------------------------------------------------------------------
M = eye(size(X,1))-X*pinv(X);

%-Get Data indices
%--------------------------------------------------------------------------
Yinds = [];
for i = 1:length(usedLabs)
    Ytmp=indchannel(S.D,usedLabs{i});
    if isempty(Ytmp)
        display([usedLabs{i} ' not found in D object'])
    end
    Yinds=[Yinds; Ytmp];
end


if (size(Yinds,1)~=size(X,1))
    error('Mismatch between data size and number of sensors with orientation information');
else
  
%-create ouput dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n');
outname = fullfile(path(S.D),['MF_' fname(S.D)]);
mfD = clone(S.D,outname);
mfD.save();

%- Work out chunk size
%--------------------------------------------------------------------------
chunkSamples= round(S.chunkSize/(8*size(S.D,1))*1e6);
begs=1:chunkSamples:size(S.D,2);
ends = (begs+chunkSamples-1);
if(ends(end)>size(S.D,2))
    ends(end)= size(S.D,2);
end
%-Run on channels needing correction
%--------------------------------------------------------------------------
vars = zeros(length(Yinds),1);

for i =1:length(begs)
    fprintf(['Modelling mean field chunk # %3.2f of %3.2f\n'],i,length(begs));
    inds = begs(i):ends(i);
    out = S.D(:,inds,:);
    Y=out(Yinds,:);
    out(Yinds,:)=M*Y;  
    mfD(:,inds,:)=out;
    vars = var(out(Yinds,:),0,2)+vars;
end

%-Update forward modelling information
%--------------------------------------------------------------------------
if (S.balance)
fprintf('Updating sensor information\n');
% TODO: Make this more fieldTrip compliant going forward.
grad = mfD.sensors('MEG');
grad.tra                = M*grad.tra;
grad.balance.previous   = grad.balance.current;
grad.balance.current    = 'mfc'; 
mfD = sensors(mfD,'MEG',grad);
% Check if any information in D.inv needs updating.
% TODO: Update to support multiple invs/forwards/modalities
if isfield(mfD,'inv')
    if isfield(mfD.inv{1},'gainmat')
        fprintf(['Clearing current forward model, please recalculate '...
            'with spm_eeg_lgainmat\n']);
        mfD.inv{1} = rmfield(mfD.inv{1},'gainmat');
    end
    if isfield(mfD.inv{1},'datareg')
        mfD.inv{1}.datareg.sensors = grad;
    end
    if isfield(mfD.inv{1},'forward')
        voltype = mfD.inv{1}.forward.voltype;
        mfD.inv{1}.forward = [];
        mfD.inv{1}.forward.voltype = voltype;
        mfD = spm_eeg_inv_forward(mfD,1);
    end
end

mfD.save();
end
%-Bad Channel Check
%--------------------------------------------------------------------------
fprintf('Checking for unusual channels\n');
SD = sqrt(vars)*1e-3;
for i = 1:length(SD)
    index= indchannel(S.D,usedLabs{i});
    if(SD(i)>S.badChanThresh)
        fprintf(['Residual on channel ' num2str(index) ', '  usedLabs{i} ': %3.2f pT\n'], SD(i));
    end
end

%-Complete
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));

end
