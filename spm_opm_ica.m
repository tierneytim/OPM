function [icD] = spm_opm_ica(S)
% Ica denoising of OPM data 
% FORMAT D = spm_opm_hfc(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.chunkSize     - max memory usage(for large datasets)           - Default 512(MB)
%   S.badChanThresh - threshold (std) to identify odd channels       - Default 2 (pT)
%   S.prefix        - prefix to filename                             - Default 'h'
%   S.reg           - regexfor selecting channels to plot            -Default 'MEG'
% Output:
%   icD               - denoised MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id: spm_opm_ica.m 8245 2022-04-26 12:59:57Z george $

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'prefix'),        S.prefix = 'ica'; end
if ~isfield(S, 'chunkSize'),     S.chunkSize = 512; end
if ~isfield(S, 'S.badChanThresh'),     S.badChanThresh = 1.5; end
if ~isfield(S, 'S.reg'),          S.reg = 'MEG'; end

%-layout creation
%--------------------------------------------------------------------------
ftfD = ftraw(S.D);
cfg=[];
[layout, ~] = ft_prepare_layout(cfg, ftfD);


regex = ['regexp_(',S.reg,')'];
displayInds = S.D.selectchannels(regex);
chanlabs = chanlabels(S.D,displayInds);
n = length(chanlabs);

chanlabs{n+1}= 'COMNT';
chanlabs{n+2}= 'SCALE';
[a,b]=spm_match_str(chanlabs',layout.label);

layout.pos=layout.pos(b,:);
layout.label=layout.label(b,:);
layout.width=layout.width(b,:);
layout.height=layout.height(b,:);


%-Run the ICA
%--------------------------------------------------------------------------

cfg                 = [];
cfg.channel         = 'all';
cfg.method          = 'fastica';
cfg.numcomponent    = 50;
cfg.feedback        = 'text';
cfg.updatesens      = 'no';
cfg.randomseed      = 454;
comp                = ft_componentanalysis(cfg, ftfD);


cfg             = [];
cfg.layout      = layout; % specify the layout file that should be used for plotting
cfg.viewmode    = 'component';
cfg.blocksize   = 10;
cfg.compscale   = 'local';
ft_databrowser(cfg, comp)


prompt="which components do you want to remove?";
comps = input(prompt);

% topo in space of grad structure 
[~,b]=spm_match_str(comp.grad.label,comp.topolabel);
topo = comp.topo(b,comps);

% topo in space of SPM structure 
[~,Yinds]=spm_match_str(comp.grad.label,chanlabels(S.D));
topo = topo(b,:);
Mgrad = eye(size(topo,1))-topo*pinv(topo);

%-Compute projector
%--------------------------------------------------------------------------
M = eye(size(topo,1))-topo*pinv(topo);

%-create ouput dataset object
%--------------------------------------------------------------------------
fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),[S.prefix fname(S.D)]);
icD = clone(S.D,outname);
icD.save();

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
trvar= zeros(length(M),size(S.D,3));
fprintf('%-40s: %30s\n','Processing Data',spm('time'));

for j=1:size(S.D,3)
    mk= S.D(Yinds,1,j);
    sk = zeros(length(Yinds),1);
    count=1;
    
    for i =1:length(begs)
        inds = begs(i):ends(i);
        
        % icD(Yinds,inds,j)=M*S.D(Yinds,inds,j) is slow (disk read)
        out = S.D(:,inds,j);
        Y=out(Yinds,:);
        out(Yinds,:)=M*Y;
        icD(:,inds,j)=out;
        
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

%-Update forward modelling information
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Updating Sensor Information',spm('time'));
grad = icD.sensors('MEG');
tmpTra= Mgrad;
%tmpTra(sinds,sinds)=M;
grad.tra                = tmpTra*grad.tra;
grad.balance.previous   = grad.balance.current;
grad.balance.current    = 'ica';
icD = sensors(icD,'MEG',grad);
% Check if any information in D.inv needs updating.
% TODO: Update to support multiple invs/forwards/modalities
if isfield(icD,'inv')
    if isfield(icD.inv{1},'gainmat')
        fprintf(['Clearing current forward model, please recalculate '...
            'with spm_eeg_lgainmat\n']);
        icD.inv{1} = rmfield(icD.inv{1},'gainmat');
    end
    if isfield(icD.inv{1},'datareg')
        icD.inv{1}.datareg.sensors = grad;
    end
    if isfield(icD.inv{1},'forward')
        voltype = icD.inv{1}.forward.voltype;
        icD.inv{1}.forward = [];
        icD.inv{1}.forward.voltype = voltype;
        icD = spm_eeg_inv_forward(icD,1);
    end
end
icD.save();

usedLabs = chanlabels(S.D,Yinds);
%-Odd Channel Check
%--------------------------------------------------------------------------
fprintf('Checking for unusual channels\n');
SD = mean(sqrt(trvar),2)*1e-3;
for i = 1:length(SD)
    index= Yinds(i);
    if(SD(i)>S.badChanThresh)
        fprintf(['Residual on channel ' num2str(index) ', '...
            usedLabs{index} ': %3.2f pT\n'], SD(i));
    end
end

%-Complete
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));
end

