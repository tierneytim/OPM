function [tHelm] = spm_affine_opreg(S)
% Read magnetometer data and optionally set up forward model
% FORMAT tHelm = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
%   S.helmetfile     - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headfile     - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.helmetref    - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headhelmetref - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.headfid
%   S.headhelmetfid
%   S.affine    
% Output:
%  tHelm       - transformed helmet object
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging

% Tim Tierney
% $Id$
spm('FnBanner', mfilename);

%-read and convert meshes to mm
%--------------------------------------------------------------------------
Native= gifti(S.headfile);

inmetres = sqrt(sum((max(Native.vertices)-min(Native.vertices)).^2))<10;
if inmetres
Native.vertices=Native.vertices*1000;
end


helmetInfo = load(S.helmetfile);
helmet = gifti(helmetInfo.Helmet_info.helmet_surface);

inmetres = sqrt(sum((max(helmet.vertices)-min(helmet.vertices)).^2))<10;
if inmetres
    helmet.vertices=helmet.vertices*1000;
end

% headhelmet = gifti('headhelmet.obj');
% inmetres = sqrt(sum((max(headhelmet.vertices)-min(headhelmet.vertices)).^2))<10;
% if inmetres
%     headhelmet.vertices=headhelmet.vertices*1000;
% end
% mesh_plot(headhelmet, helmet)


%- templates 
%--------------------------------------------------------------------------
scalp = gifti(fullfile(spm('dir'), 'canonical','scalp_2562.surf.gii'));
brain = gifti(fullfile(spm('dir'), 'canonical','cortex_5124.surf.gii'));

%- head space --> template (6 paramater)
%--------------------------------------------------------------------------
fid_template = spm_eeg_fixpnt(ft_read_headshape(fullfile(spm('dir'), 'EEGtemplates','fiducials.sfp')));
head2templatescalp = spm_eeg_inv_rigidreg(fid_template.fid.pnt(1:3,:)',S.headfid');
scalpTemplate6 = spm_mesh_transform(Native,head2templatescalp);

%- cut everything below the ears off of scanned mesh
%--------------------------------------------------------------------------
p = double(scalpTemplate6.vertices);
cut = min(fid_template.fid.pnt(1:3,3));
p(p(:,3)<cut,:) =[];

%- affine registration (12 paramater) (uniform scaling)
%--------------------------------------------------------------------------
hmm= spm_eeg_inv_icp(double(scalp.vertices'),p',[],[],[],[],S.affine);
scalpTemplate12 = spm_mesh_transform(scalpTemplate6,hmm);

%-  helmet --> head with helmet 
%--------------------------------------------------------------------------
helm2headhelm = spm_eeg_inv_rigidreg(S.headhelmetref',S.helmetref');

%-  head with helmet  --> head
%--------------------------------------------------------------------------
headhelm2head = spm_eeg_inv_rigidreg(S.headfid',S.headhelmetfid');

%-  sensors to template 
%--------------------------------------------------------------------------
sens2temp= hmm*head2templatescalp*headhelm2head*helm2headhelm;
sens=[];
sens.vertices =  helmetInfo.Helmet_info.sens_pos*1000;
senstemplate = spm_mesh_transform(sens,sens2temp);

%-  check registration
%--------------------------------------------------------------------------

p1 = [];
p1.faces = scalpTemplate12.faces;
p1.vertices = scalpTemplate12.vertices;
p2 = [];
p2.faces = scalp.faces;
p2.vertices = scalp.vertices;

f = figure()

patch(p1,'FaceColor','red','FaceAlpha',1,'EdgeColor','none')
hold on 
patch(p2,'FaceColor','blue','FaceAlpha',.3,'EdgeColor','none')

ax = gca();
camlight;
camlight(-80,-10);
hold on 
 p = senstemplate.vertices;
plot3(p(:,1),p(:,2),p(:,3),'.y','MarkerSize',20)


%-  save new sensor positions
%--------------------------------------------------------------------------
tHelm = helmetInfo.Helmet_info;
tHelm.sens_pos= senstemplate.vertices;
Helmet_info = tHelm;
[a,b,c] = fileparts(S.helmetfile);
outfile = fullfile(a,['trans_',b,c]);
save(outfile,'Helmet_info');



end