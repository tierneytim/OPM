function [vSlm] = spm_opm_vslm(S)
% solid harmonics for denoising OPM data
% FORMAT vSlm   =  spm_eeg_vslm(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.D         - SPM MEEG object      - Default: specify S.D or S.v and S.o
%   S.v         - optional positions               - Default: same as S.D
%   S.o         - optional orientations            - Default: same as S.D
%   S.or        - optional origin offset           - Default = [0,0,0]
%   S.reg       - regular or irregular (boolean)   - Default: 1
%   S.scale     - scale harmonic for stabilty      - Deult: 1
% Output:
%  vSlm            - matrix of vector spherical harmonic (n x (li^2+2*l))
%__________________________________________________________________________
%
%- example
%----------------------------------------------------------------------
% sp = spm_mesh_sphere(5);
% Nv= spm_mesh_normals(sp);
% S=[]
% S.li=2;
% S.reg=1;
% S.v=sp.vertices;
% S.o=Nv;
% H= spm_eeg_vslm(S);
%  
% 
%  
%  for i =1:size(H,2)
%  f=figure()
%  p= [];
%  p.vertices=sp.vertices;
%  p.faces=sp.faces;
%  p.EdgeColor='none';
%  col= H(:,i);
%  patch(p,'FaceVertexCData',col,'FaceColor','interp');
%  view([59,28])
%  end
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
% Tim Tierney
% $Id: spm_eeg_vslm.m 7778 2020-02-05 13:52:28Z tim $

%- handle arguments
%----------------------------------------------------------------------
if ~isfield(S, 'D')
    x=S.v(:,1);
    y=S.v(:,2);
    z=S.v(:,3);
    nx = S.o(:,1);
    ny = S.o(:,2);
    nz = S.o(:,3);
    v=S.v;
else
    p = sensors(S.D,'MEG');
    o=[p.chanori];
    
    nx = o(:,1);
    ny = o(:,2);
    nz = o(:,3);
    
    x = p.chanpos(:,1);
    y = p.chanpos(:,2);
    z = p.chanpos(:,3);
    
end

if ~isfield(S, 'or')
    S.or = [0,0,0];
end

x = (x)-S.or(1);
y = (y)-S.or(2);
z = (z)-S.or(3);
v = [x,y,z];

if ~isfield(S, 'reg')
    S.reg=1;
end
reg=S.reg;
if ~isfield(S, 'scale')
    S.scale=1;
end
%- prepare 
%----------------------------------------------------------------------

n=S.li^2+2*S.li;
r = sqrt(x.^2+y.^2+z.^2);
rbar = mean(r);
Slm = slm(v,S.li);
Slmdx=zeros(size(x,1),n);
Slmdy=zeros(size(x,1),n);
Slmdz=zeros(size(x,1),n);
vSlm=zeros(size(x,1),n);


%- Z 
%----------------------------------------------------------------------
count=1;
for l=1:S.li
    for m=-l:l
        a = sqrt((2*l+1)/4*pi* factorial(l-abs(m))/factorial(l+abs(m)));
        %a=1;
        u = m*atan2(y,x);
        um =abs(m)*atan2(y,x);
        
        if(m<0)
            Zlm=zlm(v,l,abs(m));
            t1= a*sin(um).*Zlm;
            t1(isnan(t1))=0;
            
        end
        
        if(m==0)
            Zlm=zlm(v,l,0);
            t1= a*Zlm;
            t1(isnan(t1))=0;
        end
        
        if(m>0)
            Zlm=zlm(v,l,abs(m));
            t1= a*cos(u).*Zlm;
            t1(isnan(t1))=0;
        end
        if(reg)
            Slmdz(:,count) = t1.*r.^(l)+l*z.*Slm(:,count).*r.^(l-2);
               if(S.scale)
                  Slmdz(:,count) = t1+l*z.*Slm(:,count).*r.^(-2) ;
                  Slmdz(:,count) = Slmdz(:,count).*exp(l*log(r/rbar));
                end
        else
            Slmdz(:,count) = t1./r.^(l+1)-(l+1)*z.*Slm(:,count)./r.^(l+3);
             if(S.scale)
                  Slmdz(:,count) = t1./(r.^(1))-(l+1)*z.*Slm(:,count)./(r.^(3));
                  Slmdz(:,count) = Slmdz(:,count).*exp(l*log(rbar./r));
                end
        end
        divzero = isnan(Slmdz(:,count));
        Slmdz(divzero,count)=0;
        
        count=count+1;
    end
    
end

%- y 
%----------------------------------------------------------------------
count=1;
for l=1:S.li
    for m=-l:l
        a = sqrt((2*l+1)/4*pi* factorial(l-abs(m))/factorial(l+abs(m)));
        %a=1;
        u = m*atan2(y,x);
        um =abs(m)*atan2(y,x);
        
        if(m<0)
            L = plm(z./r,l,abs(m));
            Ylm=ylm(v,l,abs(m));
            t1= a*L*abs(m).*cos(um).*x./(x.^2+y.^2);
            t1(isnan(t1))=0;
            t1= t1+a*sin(um).*Ylm;
            t1(isnan(t1))=0;     
        end
        
        if(m==0)
            Ylm=ylm(v,l,0);
            t1= a*Ylm; 
            t1(isnan(t1))=0;
        end    
            
        if(m>0)
            L = plm(z./r,l,m);
            Ylm=ylm(v,l,m);
            t1= -a*L*m.*sin(u).*x./(x.^2+y.^2);
            t1(isnan(t1))=0;
            t1= t1+a*cos(u).*Ylm;
            t1(isnan(t1))=0;
        end
          if(reg)
                Slmdy(:,count) = t1.*r.^(l)+l*y.*Slm(:,count).*r.^(l-2);
                  if(S.scale)
                  Slmdy(:,count) = t1+l*y.*Slm(:,count).*r.^(-2) ;
                  Slmdy(:,count) = Slmdy(:,count).*exp(l*log(r/rbar));
                end
          else
            Slmdy(:,count) = t1./r.^(l+1)-(l+1)*y.*Slm(:,count)./r.^(l+3);
            if(S.scale)
                  Slmdy(:,count) = t1./(r.^(1))-(l+1)*y.*Slm(:,count)./(r.^(3));
                  Slmdy(:,count) = Slmdy(:,count).*exp(l*log(rbar./r));
                end
            end  
        divzero = isnan(Slmdy(:,count));
        Slmdy(divzero,count)=0;
        count=count+1;
    end
    
end

%- X
%----------------------------------------------------------------------
count=1;
for l=1:S.li
    for m=-l:l
        a = sqrt((2*l+1)/4*pi* factorial(l-abs(m))/factorial(l+abs(m)));
        %a=1;
        u = m*atan2(y,x);
        um =abs(m)*atan2(y,x);
        if(m<0)
            L = plm(z./r,l,abs(m));
            Xlm=xlm(v,l,abs(m));
            t1= -a*L*abs(m).*cos(um).*y./(x.^2+y.^2);
            t1(isnan(t1))=0;
            t1= t1+a*sin(um).*Xlm;
            t1(isnan(t1))=0;
          
        end    
        if (m==0)
            Xlm=xlm(v,l,0);
            t1= a*Xlm;
            t1(isnan(t1))=0;
        end    
            
        if(m>0)
            L = plm(z./r,l,m);
            Xlm=xlm(v,l,m);
            t1= a*L*m.*sin(u).*y./(x.^2+y.^2);
            t1(isnan(t1))=0;
            t1=t1+a*cos(u).*Xlm;
            t1(isnan(t1))=0;
          
           
        end
          if(reg)
                Slmdx(:,count) = t1.*r.^(l)+l*x.*Slm(:,count).*r.^(l-2);
                if(S.scale)
                  Slmdx(:,count) = t1+l*x.*Slm(:,count).*r.^(-2) ;
                  Slmdx(:,count) = Slmdx(:,count).*exp(l*log(r/rbar));
                end
            else
                Slmdx(:,count) = t1./(r.^(l+1))-(l+1)*x.*Slm(:,count)./(r.^(l+3));
                if(S.scale)
                  Slmdx(:,count) = t1./(r.^(1))-(l+1)*x.*Slm(:,count)./(r.^(3));
                  Slmdx(:,count) = Slmdx(:,count).*exp(l*log(rbar./r));
                end
            end
        divzero = isnan(Slmdx(:,count));
        Slmdx(divzero,count)=0;
        
       count=count+1;
    end
    
end


%- cleanup
%----------------------------------------------------------------------
for i = 1:n
vSlm(:,i)=Slmdz(:,i).*nz+Slmdy(:,i).*ny+Slmdx(:,i).*nx;
end

end

function [Xlm] = xlm(v,l,m)
x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);
b= (-1)^m * 2^l;
Xlm=0;

for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    num = -x.*z.^(k-m).*(k-m).*(x.^2+y.^2).^(m/2) + (m*x.*z.^(k-m+2)).*(x.^2+y.^2).^((m-2)/2);
    %numa= (x.*z.^(k-m)).*((x.^2+y.^2).^(m/2-1));
    %numb=m*x.^2+m*y.^2+m*z.^2-k*x.^2-k*y.^2;
    %Xlm=Xlm+b*c*(numa.*numb)./(r.^(2+k));
    Xlm=Xlm+b*c*(num)./(r.^(2+k));
    
end

end

function [Ylm] = ylm(v,l,m)
x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);
b= (-1)^m * 2^l;
Ylm=0;

for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    num = -y.*z.^(k-m).*(k-m).*(x.^2+y.^2).^(m/2) + (m.*y.*z.^(k-m+2)).*(x.^2+y.^2).^((m-2)/2);
    %numa = (y.*z.^(k-m)).*((x.^2+y.^2).^(m/2-1));
    %numb =m*x.^2+m*y.^2+m*z.^2-k*x.^2-k*y.^2;
   % Ylm=Ylm+b*c*(numa.*numb)./(r.^(2+k));
    Ylm=Ylm+b*c*(num)./(r.^(2+k));
end

end

function [Zlm] = zlm(v,l,m)
x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);
b= (-1)^m * 2^l;
Zlm=0;

for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    num = z.^(k-m-1).*(k-m).*(x.^2+y.^2).^((m+2)/2) + (-m*z.^(k-m+1)).*(x.^2+y.^2).^(m/2);
    %numa = (z.^(k-m-1)).*((x.^2+y.^2).^(m/2));
    %numb =m*x.^2+m*y.^2+m*z.^2-k*x.^2-k*y.^2;
    %Zlm=Zlm+b*c*(-numa.*numb)./(r.^(2+k));
    Zlm=Zlm+b*c*(num)./(r.^(2+k));
    
    end

end

function [Slm] = slm(v,li)

x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);

n=li^2+2*li;
Slm=zeros(size(v,1),n);
count=1;

for l=1:li
    for  m = -l:l
        a = sqrt((2*l+1)/4*pi* factorial(l-abs(m))/factorial(l+abs(m)));
       % a=1;
        if(m<0)
            L = plm(z./r,l,abs(m));
            %Slm(:,count) =  a*L.*sin(abs(m)*acos(x./sqrt(x.^2+y.^2)));
            Slm(:,count) =  a*L.*sin(abs(m)*atan2(y,x));
        elseif m==0
            L = plm(z./r,l,0);
            Slm(:,count) =a*L;
        else
            L = plm(z./r,l,m);
            %Slm(:,count) =  a*L.*cos(m*acos(x./sqrt(x.^2+y.^2)));
            Slm(:,count) =  a*L.*cos(m*atan2(y,x));
        end
        divzero = isnan(Slm(:,count));
        Slm(divzero,count)=0;
        count=count+1;
    end
end

end

function [pl] = plm(x,l,m)
b= (-1)^m * 2^l;
pl =0;
for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    pl=pl+(b*(1-x.^2).^(m/2))*c.*x.^(k-m);
end
end
