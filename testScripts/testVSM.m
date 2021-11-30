% sphere to evaluate harmonics on
sp = spm_mesh_sphere(5);

% general solution
S=[];
S.li=3;
S.reg=1;
S.v=sp.vertices*50;
S.o=sp.vertices;
H2= spm_opm_vslm(S);

% radial orienation and position ar same on unit sphere
ex = sp.vertices(:,1);
ey = sp.vertices(:,2);
ez = sp.vertices(:,3);

% position on scaled sphere
x = sp.vertices(:,1)*50;
y = sp.vertices(:,2)*50;
z = sp.vertices(:,3)*50;


% brute force solution 
H=[];
L1 = [ey, ez,ex ];
H= [H,L1];
L2 =[y.*ex+x.*ey, z.*ey+y.*ez, ...
    -x.*ex-y.*ey+2*z.*ez, ...
    z.*ex+x.*ez, ...
    x.*ex-y.*ey];
H= [H,L2];
L3 = [(2.*x.*y).*ex+(x.^2-y.^2).*ey,...
    (y.*z).*ex+(x.*z).*ey + (x.*y).*ez,...
    (2*x.*y).*ex+(x.^2+3*y.^2-4*z.^2).*ey+(-8*z.*y).*ez,...
    (-6*x.*z).*ex+(-6*y.*z).*ey+(6*z.^2-3*x.^2-3*y.^2).*ez,...
    (3*x.^2+y.^2-4*z.^2).*ex+(2*x.*y).*ey+(-8*x.*z).*ez,...
    (2*x.*z).*ex+(-2*y.*z).*ey+(x.^2-y.^2).*ez,...
    (3*y.^2-3*x.^2).*ex+(6*x.*y).*ey];
H= [H,L3];

% correlate general and brute force solution 
 cs = zeros(size(H2,2),1);
 for i=1:size(H2,2)
     cs(i)=abs(corr(H2(:,i),H(:,i)));
 end
 
 % if both are perfectly corelated sum is 15 (15 regresors)
 if (sum(cs)==15)
     disp('VSM Test Passed');
 else
     disp('VSM Test Failed');
 end