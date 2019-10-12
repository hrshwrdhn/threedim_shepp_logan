function [p]=phantom3kak(nx,ny,nz)
% use to generate three dimension phantom

%      
%         A      a     b     c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e =    [  2  .6900  .920  .900      0       0       0      0      0      0
        -.98  .6624  .874  .880      0       0       0      0      0      0
        -.02  .4100  .160  .210   -.22       0    -.25    108      0      0
        -.02  .3100  .110  .220    .22       0    -.25     72      0      0
         0.02  .0460  .046  .046      0     .1    -.25      0      0      0
         .02  .0460  .046  .046      0      .1    -.25      0      0      0
         .01  .0460  .023  .020   -.8    -.65    -.25      0      0      0
         .01  .0460  .023  .020    .06    -.065    -.25     90      0      0
         .02  .560  .040  .100    .06   -.105    .625     90      0      0
        -.02  .0560  .056  .100      0    .100    -.625      0      0      0 ];



n=128;
nx=n;
ny=n;
nz=n;
p = zeros([nx ny nz]);
rngx= ( (0:nx-1)-(nx-1)/2 ) / ((nx-1)/2); 
rngy = -1*( (0:ny-1)-(ny-1)/2 ) / ((ny-1)/2); 
rngz = -1*( (0:nz-1)-(nz-1)/2 ) / ((nz-1)/2);

[x,y,z] = meshgrid(rngx,rngy,rngz);
coord = [flatten(x); flatten(y); flatten(z)];
p = flatten(p);

for k = 1:10    
   A = e(k,1);            % Amplitude change for this ellipsoid
   asq = e(k,2)^2;        % a^2
   bsq = e(k,3)^2;        % b^2
   csq = e(k,4)^2;        % c^2
   x0 = e(k,5);           % x offset
   y0 = e(k,6);           % y offset
   z0 = e(k,7);           % z offset
   phi = e(k,8)*pi/180;   % first Euler angle in radians
   theta = e(k,9)*pi/180; % second Euler angle in radians
   psi = e(k,10)*pi/180;  % third Euler angle in radians
   
   cphi = cos(phi);
   sphi = sin(phi);
   ctheta = cos(theta);
   stheta = sin(theta);
   cpsi = cos(psi);
   spsi = sin(psi);
   
   % Euler rotation matrix
   alpha = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
            -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
            stheta*sphi                  -stheta*cphi                ctheta];        
   
   % rotated ellipsoid coordinates
   coordp = alpha*coord;
   
   idx = find((coordp(1,:)-x0).^2./asq + (coordp(2,:)-y0).^2./bsq + (coordp(3,:)-z0).^2./csq <= 1);
   p(idx) = p(idx) + A;

end

p = reshape(p,[n n n]);




return;

function out = flatten(in)
out = reshape(in,[1 prod(size(in))]);
return;
   
