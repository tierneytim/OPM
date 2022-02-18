function [T] = spm_opm_sourcepow(S)
% computed a sliding window of t-stats(power changes from baseline)in source space  
% FORMAT D = spm_opm_sourcepow(S)
%   S               - input structure
%  fields of S:
%   S.D           - epoched, filtred SPM MEEG object   - Default: no Default
%   S.basewin     - size of window in ms eg [-1500 0]  - Default: no Default              
%   S.actwin      - size of window in ms eg [0 1500]   - Default: no Default              
% Output:
%   pD               - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging
%
%
% Tim Tierney
% $Id$

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),      error('an MEEG object must be supplied.'); end
if ~isfield(S, 'basewin'),    error('a time window(ms) must be supplied.'); end
if ~isfield(S, 'actwin'),    error('a time window(ms) must be supplied.'); end

%- Get source model
%--------------------------------------------------------------------------
nverts = length(S.D.inv{1}.forward.mesh.vert);
[~,S.D] = spm_eeg_lgainmat(S.D,1:nverts);
S.D.save();

lgain = load(S.D.inv{1}.gainmat);
indChans = indchannel(S.D,lgain.label);
L=lgain.G;
brain=gifti(S.D.inv{1}.mesh.tess_ctx);

%- Get inversion(min norm here but does not necessarily have to be )
%--------------------------------------------------------------------------
Nsens=size(L,1);
SNR2=100;
lambda2=trace(L*L')/(Nsens*SNR2);
temp=L*L'+lambda2*eye(Nsens);
Linv=L'/temp;


%- time window baseline
%--------------------------------------------------------------------------
baseinds = find(S.D.time()>=S.basewin(1,1)/1000 & S.D.time<=S.basewin(1,2)/1000);
brainbase = zeros(size(Linv,1),size(S.D,3));
for i = 1:size(S.D,3)
    cb = cov(S.D(indChans,baseinds,i)');
    for j = 1:size(brainbase,1)
        li = Linv(j,:);
        lt= li';
        brainbase(j,i) = li*cb*lt;       
    end
end
brainbase2= (brainbase).^(1/4);

%- time window activation
%--------------------------------------------------------------------------
T = zeros(size(L,2),size(S.actwin,1));
for k = 1:size(S.actwin,1)
    k
    interestinds = find(S.D.time()>=S.actwin(k,1)/1000 & S.D.time<=S.actwin(k,2)/1000);
    braininterest = zeros(size(Linv,1),size(S.D,3));
    for i = 1:size(S.D,3)
        ci = cov(S.D(indChans,interestinds,i)');
        for j = 1:size(brainbase,1)
            li = Linv(j,:);
            lt= li';
            braininterest(j,i) = li*ci*lt;
        end
    end
    
    
    braininterest2= (braininterest).^(1/4);
    delta = braininterest2-brainbase2;
    
    Z = mean(delta,2);
    se = std(delta,[],2)/sqrt(size(S.D,3));
    T(:,k)=  Z./se;
end

maT = max(max(abs(T)));
Tplot = T;
Tplot(4336,:)=maT;
Tplot(4686,:)=-maT;


%for i = 1:size(T,2)
ax= spm_mesh_render(brain);
ax = spm_mesh_render('Overlay',ax,Tplot(:,1));
ax = spm_mesh_render('View',ax, [-279 18]);
ax = spm_mesh_render('ColourBar',ax,'on');
ax = spm_mesh_render('ColourMap',ax,colormap(flipud(brewermap(64,'RdBu'))));
ax.colourbar.Subtitle.String = 't-Statistic';
ax.colourbar.Subtitle.FontSize=16;
b = uicontrol('Parent',ax.figure,'Style','slider','Position',[81,54,419,23],...
              'value',1, 'min',1, 'max',size(T,2));
bl3 = uicontrol('Parent',ax.figure,'Style','text','Position',[240,25,100,23],...
                'String','Time','BackgroundColor',ax.figure.Color);
            
b.Callback = @(es,ed)spm_mesh_render('Overlay',ax,Tplot(:,round(es.Value))); 

end




