function outs = poiGridSurface(surface,space,lowerThresh,upperThresh)
% Find where a grid intersects a surface (ideally convex)
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

% Create Grid 
%--------------------------------------------------------------------------
bbMin = min(surface.vertices);
bbMax = max(surface.vertices);
x = (bbMin(1)-space):space:(bbMax(1)+space);
y = (bbMin(2)-space):space:(bbMax(2)+space);
z = (bbMin(3)-space):space:(bbMax(3)+space);
[X,Y,Z] = ndgrid(x,y,z);
grid   = [X(:) Y(:) Z(:)];

% Keep only outer layer of grid
%--------------------------------------------------------------------------
mi = min(grid);
ma = max(grid);

miX = grid(find(abs(grid(:,1) - mi(1))<=.01),:);
miY = grid(find(abs(grid(:,2) - mi(2))<=.01),:);
miZ = grid(find(abs(grid(:,3) - mi(3))<=.01),:);

maX = grid(find(abs(grid(:,1) - ma(1))<=.01),:);
maY = grid(find(abs(grid(:,2) - ma(2))<=.01),:);
maZ = grid(find(abs(grid(:,3) - ma(3))<=.01),:);




% ny
%--------------------------------------------------------------------------
outs = zeros(length(maY),3);

for j = 1:length(maY)
Ia = maY(j,:);
Ib = maY(j,:)+200.*[0,-1,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outsny = outs(ind,:);
display('Negative Y')

% py
%--------------------------------------------------------------------------
outs = zeros(length(miY),3);

for j = 1:length(miY)
Ia = miY(j,:);
Ib = miY(j,:)+200.*[0,1,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outspy = outs(ind,:);
display('Positive Y')

% nx
%--------------------------------------------------------------------------
outs = zeros(length(maX),3);

for j = 1:length(maX)
Ia = maX(j,:);
Ib = maX(j,:)+200.*[-1,0,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outsnx = outs(ind,:);
display('Negative X');

% px
%--------------------------------------------------------------------------
outs = zeros(length(miX),3);

for j = 1:length(miX)
Ia = miX(j,:);
Ib = miX(j,:)+200.*[1,0,0];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outspx = outs(ind,:);
display('Postive X')

% nz
%--------------------------------------------------------------------------
outs = zeros(length(maZ),3);

for j = 1:length(maZ)
Ia = maZ(j,:);
Ib = maZ(j,:)+200.*[0,0,-1];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outsnz = outs(ind,:);

display('Negative Z');

% pz
%--------------------------------------------------------------------------
outs = zeros(length(miZ),3);

for j = 1:length(miZ)
Ia = miZ(j,:);
Ib = miZ(j,:)+200.*[0,0,1];

m = zeros(3,3);
m(:,1) = Ia-Ib;

nfaces = length(surface.faces);
tuv = zeros(nfaces,3);
verts = double(zeros(3,3));
Y = zeros(3,1);

for i = 1:nfaces
verts = surface.vertices(surface.faces(i,:),:);
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
outspz = outs(ind,:);

display('Positive Z');
% Initial Points
%--------------------------------------------------------------------------
uspz = unique(round(outsnz,3),'rows');
% nn= [];
% rem = zeros(length(uspz),1);
% dist = squareform(pdist(uspz));
% 
% for i =1:length(uspz)
%     sdist= sort(dist(:,i));
%     nn= min(sdist(2:4));
%     if(nn<upperThresh)
%         rem(i)=1;
%     else
%         rem(i)=0;
%     end
%         
% end
% uspz = uspz(boolean(rem),:);

% Constrained adding of points
%--------------------------------------------------------------------------
other = [unique(round(outsnx,3),'rows');...
         unique(round(outspy,3),'rows');...
         unique(round(outsny,3),'rows');...
         unique(round(outspz,3),'rows');...
         unique(round(outspx,3),'rows');...
         ];
    us = uspz;
    
for i = 1:length(other)
ot = other(i,:);
di = sqrt(sum((bsxfun(@minus,us,ot)).^2,2));
if(all(di>lowerThresh))
   us = [us;ot];
   
% rem = zeros(length(us),1);
% dist = squareform(pdist(us));
% for i =1:length(us)
%     sdist= sort(dist(:,i));
%     nn= min(sdist(2:6));
%     if(nn>upperThresh)
%         rem(i)=0;
%     else
%         rem(i)=1;
%     end
%         
% end 
end

end

% return
%--------------------------------------------------------------------------     
outs = us;

end
