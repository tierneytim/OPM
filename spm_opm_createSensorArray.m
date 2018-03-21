function [pos,ori] = spm_opm_createSensorArray(S)
% Given a scalp surface even samples the surface with sensors 
% 
% FORMAT [pos,ori] = spm_opm_createSensorArray(S)
%   S               - input structure
% Fields of S:
%   S.D            - SPM M/EEG object with surface meshes 
%   S.offset       - distance to place sensors(mm) from scalp surface
%   S.space        - distance between sensors(mm) 
%   S.wholehead    - boolean: Should whole scalp surface should be covered?
% _________________________________________________________________________
% Copyright (C) <2017>  <Tim Tierney>
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%__________________________________________________________________________

% Args
%--------------------------------------------------------------------------
D = S.D;
offset = S.offset;       
space = S.space;
wholehead = S.wholehead;
% Meshes
%--------------------------------------------------------------------------
scalp = gifti(D.inv{1}.mesh.tess_scalp);
cortex = gifti(D.inv{1}.mesh.tess_ctx);
Nf = size(scalp.faces,1);
lp = min(cortex.vertices(:,3));

% Face Positions(Downsampled)
%--------------------------------------------------------------------------

if(~wholehead)
    C= scalp.vertices(:,3)>(lp);
dscalp = spm_mesh_split(scalp,C);
else 
    dscalp =scalp;
end
p = [];
p.faces =dscalp.faces;
p.vertices = double(dscalp.vertices);
dscalp = reducepatch(p,.05);

[Nv,Nf]= spm_mesh_normals(dscalp,true);
posi = zeros(size(dscalp.faces));
cog = mean(scalp.vertices);

for i = 1:length(posi)
    whichVerts = dscalp.faces(i,:);
    cos = dscalp.vertices(whichVerts,:);
    posi(i,:) = mean(cos);   
end

% Face Orientations (Downsampled)
%--------------------------------------------------------------------------
ori = zeros(size(Nf));

for i = 1:length(Nf)
ad = posi(i,:) + 5*Nf(i,:);
subtrac = posi(i,:) - 5*Nf(i,:);

d1 = sum((cog - ad).^2);
d2 = sum((cog - subtrac).^2);

if(d2>d1)
    ori(i,:) = -Nf(i,:);
else
    ori(i,:) = Nf(i,:);
end
end

% Evenly sampled Positions on Downsampled scalp
%--------------------------------------------------------------------------
outs = poiGridSurface(dscalp,space,space*.8,space*1.25);

% Find what face we're on
%--------------------------------------------------------------------------
faceInd = zeros(size(outs,1),1);
for i = 1:length(outs)
    [mi,ind] = min(sqrt(sum(bsxfun(@minus,posi,outs(i,:)).^2,2)));
    faceInd(i)=ind;

end

% Project to used scalp
%--------------------------------------------------------------------------
pos = poiLineSurface(scalp,outs,ori(faceInd,:),space);


% Check if wholehead is requested
%--------------------------------------------------------------------------
if(~wholehead)
    C= scalp.vertices(:,3)>(lp);
    scalp = spm_mesh_split(scalp,C);
end
% Face Positions
%--------------------------------------------------------------------------
[Nv,Nf]= spm_mesh_normals(scalp,true);
posi = zeros(size(scalp.faces));
cog = mean(scalp.vertices);

for i = 1:length(posi)
    whichVerts = scalp.faces(i,:);
    cos = scalp.vertices(whichVerts,:);
    posi(i,:) = mean(cos);   
end

% Face Orientations 
%--------------------------------------------------------------------------
ori = zeros(size(Nf));

for i = 1:length(Nf)
ad = posi(i,:) + 5*Nf(i,:);
subtrac = posi(i,:) - 5*Nf(i,:);

d1 = sum((cog - ad).^2);
d2 = sum((cog - subtrac).^2);

if(d2>d1)
    ori(i,:) = -Nf(i,:);
else
    ori(i,:) = Nf(i,:);
end
end


% Find what face we're on
%--------------------------------------------------------------------------
faceInd = zeros(size(pos,1),1);
for i = 1:length(pos)
    [mi,ind] = min(sqrt(sum(bsxfun(@minus,posi,pos(i,:)).^2,2)));
    faceInd(i)=ind;
end


% add points if faces are empty and not near a sensor
%--------------------------------------------------------------------------

for i = 1:length(posi)
p=posi(i,:);
di= sort(sqrt(sum(bsxfun(@minus,pos,p).^2,2)));
if(di(1)>(space))
    pos = [pos;p];
    faceInd= [faceInd;i];
end
end
% return
%--------------------------------------------------------------------------

ori = ori(faceInd,:);
pos = pos+ori*offset;
end