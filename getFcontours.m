function Uyx = getFcontours(Uyx,Fyx,th,k)
% function Uyx = plotcontours(Uyx,Fyx,y,x,th,k)
%   th = threshold
%   k = k^th color choice

rgb = [
    255 255 255
    0 255 0
    255 255 0
    255 0 0
    255 0 255    %YL: megenta for the specified contour
    ];

dth = th/5;
th1 = th-dth;
th2 = th+dth;

[Ny Nx Nc] = size(Uyx);
uy = 2:Ny-1;
ux = 2:Nx-1;

mask = double(Fyx(uy,ux,1)>th1 & Fyx(uy,ux,1)<th2);
Uyx(uy,ux,1) = uint8(mask*rgb(k,1)) + uint8(~mask).*Uyx(uy,ux,1);  % red
Uyx(uy,ux,2) = uint8(mask*rgb(k,2)) + uint8(~mask).*Uyx(uy,ux,2);   % green
Uyx(uy,ux,3) = uint8(mask*rgb(k,3)) + uint8(~mask).*Uyx(uy,ux,3);   % blue



