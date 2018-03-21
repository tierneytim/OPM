function [D] = spm_opm_removeOutlierTrials(S)
% Removes trials that are marked as outlers
% FORMAT D = spm_opm_removeOutlierTrials(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object              - Default: no Default
%   S.thresh        - No. of MADs for outlier      - Default: 3
%                      declaration                 
%   S.prefix        - File prefix                  - Default: o 
%   S.exclude       - nx1 cell array of channels   - Default: []
%                       to be excluded from 
%                       variance calcualtions
% Output:
%  D           -  MEEG object (also written to disk)
%
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging


% set Defaults
%--------------------------------------------------------------------------
if ~isfield(S, 'D'),        error('D is required'); end
if ~isfield(S, 'thresh'),   S.thresh=3;end
if ~isfield(S, 'prefix'),   S.prefix='o';end
if ~isfield(S, 'exclude'),  S.exclude={};end

% Exclude channels if necessary
% -------------------------------------------------------------------------
indExclude = selectchannels(S.D,S.exclude);
megInd = S.D.indchantype('MEG');
megInd = setdiff(megInd,indExclude);
megChans = S.D(megInd,:,:);

% Trial by Trial variance estimate accross channels
% -------------------------------------------------------------------------
varTime= squeeze(std(megChans,0,2));
meanTrialVar= mean(varTime,1);

% Medians and MADs
%--------------------------------------------------------------------------
med= median(meanTrialVar);
ma= mad(meanTrialVar,1);
thresh = med+ma*1.4826*S.thresh;

% Threshold and find
%-------------------------------------------------------------------------- 
badBin = meanTrialVar>thresh;
nBad= sum(badBin);
retain = find(~badBin);

% Crop object
%--------------------------------------------------------------------------
[a,b,c]= fileparts(fullfile(S.D));
fname = fullfile(a,[S.prefix,b,c]);
args=[];
args.D=S.D;
args.inds = retain;
args.fname= fname;
D = spm_opm_selectTrials(args);
% summary Figure
%--------------------------------------------------------------------------
figure();
stem(meanTrialVar);
hold on
plot(find(badBin),meanTrialVar(badBin),'.r','MarkerSize', 30);
lthresh = repmat(thresh,length(meanTrialVar),1);
plot(1:length(meanTrialVar),lthresh,'--r','LineWidth',3)
ylabel('Trial Standard Deviation (fT)')
xlabel('Trial No.')
set(gcf,'Color','w')
legend('Standard Deviation','Outliers','Threshold')

end