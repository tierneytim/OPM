function [SSPD,Yinds] = spm_opm_SSP(S)
% Estimates 3 eigenmodes from training data and then removes from test data(SSP)
% FORMAT D = spm_opm_SSP(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.train         - SPM MEEG training  object                      - Default: uses first T seconds of test data for training
%   S.T             - number of seconds to compute eigenmodes from   - Defualt: 5
%   S.chunkSize     - max memory usage(for large datasets)           - Default 512(MB)
% Output:
%   SSPD            - denoised MEEG object (also written to disk)
%   Yinds           - the indices of filtered channels
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
%
% The code here implements a very basic Signal Space Projection. No forward
% model corrections are currently made and therefore this function should 
% not to be used for general denoising of  brain data. Function is useful
% for quickly assessing white nose floor of OPM arrays when sensor orientations 
% are not available and spm_opm_mfc cannot be used. In all other cases 
% spm_opm_mfc shuold be  used 

% Tim Tierney
% $Id: spm_opm_SSP.m 7646 2019-07-25 13:58:46Z tim $

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'train'),         S.train=S.D; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'balance'),       S.balance = 0; end
if ~isfield(S, 'usebadchans'),   S.usebadchans = 0; end
if ~isfield(S, 'T'),             S.T = 5; end
if ~isfield(S, 'badChanThresh'), S.badChanThresh = 50; end

%-Get design matrix
%--------------------------------------------------------------------------
n=3;
Yinds = indchantype(S.train,'MEGMAG')';
bcInds = badchannels(S.D);
if(S.usebadchans)
    Yinds = Yinds;    
else
    Yinds = setdiff(Yinds,bcInds);
end
usedLabs = chanlabels(S.train,Yinds);
testY=S.D(Yinds,find(S.train.time()<S.T),:);
G= testY*testY';
[U,~,~]= svd(G);
X=U(:,1:n);
fprintf('%-40s: %30s\n','Created Design Matrix',spm('time'));

%-Compute projector
%--------------------------------------------------------------------------
M = eye(size(X,1))-X*X';

%-create ouput dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n');
outname = fullfile(path(S.D),['SSP_' fname(S.D)]);


SSPD = clone(S.D,outname);
SSPD.save();

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
trvar= zeros(length(Yinds),size(S.D,3));
fprintf('%-40s: %30s\n','Processing Data',spm('time'));

for j=1:size(S.D,3)
    mk= S.D(Yinds,1,j);
    sk = zeros(length(Yinds),1);
    count=1;
    
    for i =1:length(begs)
        inds = begs(i):ends(i);
        
        % mfD(Yinds,inds,j)=M*S.D(Yinds,inds,j) is slow (disk read)
        out = S.D(:,inds,j);
        Y=out(Yinds,:);
        out(Yinds,:)=M*Y;
        SSPD(:,inds,j)=out;
        
        % accurate running variance for identifying odd channels
        % (https://www.johndcook.com/blog/standard_deviation/)
        for l = 1:length(inds)
            xk = out(Yinds,l);
            mkprev = mk;
            mk = mkprev +(xk-mkprev)/count;
            sk=sk+(xk-mkprev).*(xk-mk) ;
            count=count+1;
        end
        
    end
    trvar(:,j)=sk/(count-1);
end
SSPD.save();


%-Bad Channel Check
%--------------------------------------------------------------------------
fprintf('Checking for unusual channels\n');
SD = mean(sqrt(trvar),2)*1e-3;
for i = 1:length(SD)
    index= indchannel(S.D,usedLabs{i});
    if(SD(i)>S.badChanThresh)
        fprintf(['Residual on channel ' num2str(index) ', '...
            usedLabs{i} ': %3.2f pT\n'], SD(i));
    end
end

%-Complete
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));

end
