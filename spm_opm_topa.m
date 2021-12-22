function [M] = spm_opm_topa(S)
% Reconstruct sensor data on scalp surface
% FORMAT D = spm_opm_topa(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                                - Default: no Default
%   S.orientation   - 'RAD' or 'TAN'                                 - Default: 'TAN'
%   S.FWHM          - smoothing(mm) to reconstruct on scalp          - Default: 70
%   S.plot          - flag to plot topography                        - Deafult: 1
%   S.timepoint     - timepoint to plot(ms)                          - Deafult: no Default
%   S.T             - values for each sensor                         - Deafult: []
% Output:
%   M               - Patch object
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id$


%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),             error(errorMsg); end
if ~isfield(S, 'orientation'),   S.orientation = 'TAN'; end
if ~isfield(S, 'FWHM'),          S.FWHM = 70; end
if ~isfield(S, 'plot'),          S.plot = 1; end
errorMsg = 'A timepoint needs to be supplied if S.T is empty';
if ~isfield(S, 'timepoint') && ~isfield(S, 'T'), error(errorMsg); end

%-create matrix to hold positions and values
%--------------------------------------------------------------------------
s= sensors(S.D,'MEG');
labs= s.label;
Inds = contains(labs,S.orientation);
Labs = labs(Inds);
Mat =zeros(sum(Inds),4);
Yinds = indchannel(S.D,Labs);

%-Check if specific timepoint requested
%--------------------------------------------------------------------------
try
    t= S.D.time();
    [~,tind]=min(abs(t*1000-S.timepoint));
    T = S.D(Yinds,tind);
catch
    T=S.T;
end
Mat(:,1) = T;
Mat(:,2:4)=s.chanpos(Inds,:);

%- process scalp and work out how many smoothing iterations are required
%--------------------------------------------------------------------------
scalp = gifti(S.D.inv{1}.mesh.tess_scalp);
v=scalp.vertices;
[~,Di] = spm_mesh_neighbours(scalp,1);
muNeighbour= mean(mean(Di));
n= round((S.FWHM/muNeighbour)^2);

%- place supplied values at closest scalp location and smooth
%--------------------------------------------------------------------------

col =zeros(length(v),1);

for i =1:size(Mat,1)
    [~,ind]=min(sum(bsxfun(@minus,v,Mat(i,2:4)).^2,2));
    col(ind)=Mat(i,1);
end
col =spm_mesh_smooth(scalp,col,n);

%- create patch object
%--------------------------------------------------------------------------
scale=max(Mat(:,1))/max(col);

M= [];
M.vertices=scalp.vertices;
M.faces=scalp.faces;
M.EdgeColor='none';
M.FaceVertexCdata=col.*scale;
M.FaceColor='interp';

if S.plot
    figure()
    M=patch(M);
    colormap('jet')
    fig= gcf;
    fig.Color=[1,1,1];
    c=colorbar;
    c.Label.String= 'Magnetic Flux Density (fT)';
    ax=gca();
    ax.Visible='off';
end


end
