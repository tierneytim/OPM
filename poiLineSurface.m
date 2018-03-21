function outs = poiLineSurface(surface,points,normals,dist)
% Find where a point intersects a surface 
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
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging


%Point of intersection
%--------------------------------------------------------------------------
outs = zeros(length(points),3);

% compute all face centres
posi = zeros(size(surface.faces));
for i = 1:length(posi)
    whichVerts = surface.faces(i,:);
    cos = surface.vertices(whichVerts,:);
    posi(i,:) = mean(cos);   
end


for j = 1:length(points)
Ia = points(j,:);
Ib = points(j,:)+dist.*normals(j,:);

m = zeros(3,3);
m(:,1) = Ia-Ib;

di = sqrt(sum(bsxfun(@minus,posi,Ia).^2,2));
faceInd = find(di<dist);


smallFaces = surface.faces(faceInd,:);
nfaces = length(smallFaces);

tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(smallFaces(i,:),:);
m(:,2) = verts(2,:)-verts(1,:);
m(:,3) = verts(3,:)-verts(1,:);
Y = (Ia-verts(1,:))';
tuv(i,:) = m\Y;
end

a = tuv(:,1)< 1 & tuv(:,1)>0;
b = tuv(:,2)< 1 & tuv(:,2)>0;
c = tuv(:,3)< 1 & tuv(:,3)>0;
d = (tuv(:,2)+tuv(:,3))<=1;
e = a&b&c&d;
mul =min(tuv(e,1));

if(~isempty(mul))
outs(j,:) = Ia+(Ib-Ia)*mul;
end

end
ind = sum(abs(outs),2)>0;
outs = outs(ind,:);
display('Positive Y')



end
