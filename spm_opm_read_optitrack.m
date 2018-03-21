function [opti] = spm_opm_read_optitrack(S)
% FORMAT [opti] = spm_opm_read_lvm(S)
%   S               - input structure
% Optional fields of S:
%   S.filename      - filepath to csv file             -Default: no Default                      
%   S.headerlength  - integer specifying how many      -Default: 23
%                     lines of file are header
%   S.avPos         - flag to average pos              -Deafult: 1
%   S.optiTrigger   - optiTrack Trigger                -Deafult: no default

% Output: opto - output Structure
%  Fields of opto:
%   opti.motion     - motion data
% _________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'fname'),        error('filename needs to be provided');end
if ~isfield(S, 'headerlength'), S.headerlength = 8; end
if ~isfield(S, 'avPos'),        S.avPos = 1; end

%-find start of data
%--------------------------------------------------------------------------

data = dlmread(S.fname, ',',8,0);
fs = 1/mean(diff(data(:,2)));
rot=data(:,3:5);
posi=data(:,6:end);

if(S.avPos)
    x= mean(posi(:,1:3:31),2);
    y= mean(posi(:,2:3:32),2);
    z= mean(posi(:,3:3:33),2);
    pos=[x,y,z];
else
    pos=posi;
end

%- Interpolate to other space
%--------------------------------------------------------------------------

%-Output Struct
%--------------------------------------------------------------------------
opti.motion= [pos,rot];
opti.fs= fs;
opti.time=data(:,2);
end