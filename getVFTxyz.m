
function [VVOL TT FF FFF] = getVFTxyz(F,VOL,T,conFLAG,conVALUE,ix,iy,iz,Nx,Ny,Nz,dx,dy,dz,wt)

%%YL: check on the lookOGSoutput.m for the meaning of the input arguments 
%
% % view indices
% ix = round(xs/dx + Nx/2+1/2);
% iy = round(ys/dy + Ny/2+1/2);
% iz = round(zs/dz + 1);
% ixs = ix; iys = iy; izs = iz; % photon launch point
% % x,y,z
% x = ([1:Nx]-Nx/2-1/2)*dx;
% y = ([1:Ny]-Ny/2-1/2)*dx;
% z = ([1:Nz]-1)*dz;

%%% Tissue Image (RGB image)
Vzx = squeeze(VOL(iy,:,:))'; % in z,x plane through source
Vzy = squeeze(VOL(:,ix,:))';
Vyx = squeeze(VOL(:,:,iz)); % in z,x plane through source
maxV = max(VOL(:));
Vzx = Vzx/maxV*255;
Vzy = Vzy/maxV*255;
Vyx = Vyx/maxV*255;
for j=1:3  % create RGB images
    Uzy(:,:,j) = uint8(Vzy);
    Uzx(:,:,j) = uint8(Vzx);
    Uyx(:,:,j) = uint8(Vyx);
end

% composite RGB image
Nyy  = 2*wt+Nz+Ny+50;   % total height of composite image
Nxx  = 2*wt+2*Nx;      % total  width of composite image

%%% gray image
VVOL = zeros(Nyy,Nxx,3,'uint8');
for j=1:3, VVOL(:,:,j) = uint8(100);end % gray background
VVOL(1:Nz,wt+(1:Ny),:)            	= (Uzy);
VVOL(wt+Nz +(1:Nx),wt+(1:Nx),:)     = (Uzx);
VVOL(wt+Nx+(1:Ny),2*wt+Nx+(1:Nx),:)	= (Uyx);

%%% iso-Fluence contours
Fzx = squeeze(F(iy,:,:))'; % in z,x plane through source
Fzy = squeeze(F(:,ix,:))';
Fyx = squeeze(F(:,:,iz)); % in z,x plane through source
FF  = zeros(Nyy,Nxx,3,'uint8');

FFF = zeros(Nyy,Nxx);
FFF(1:Nz,wt+(1:Ny))              	= (Fzy);
FFF(wt+Nz+(1:Nx),wt+(1:Nx))      	= (Fzx);
FFF(wt+Nx+(1:Ny),2*wt+Nx+(1:Nx))	= (Fyx);

for j=1:3, FF(:,:,j) = uint8(100);end % gray background
if conFLAG == 0
    tth = [0.1 1 10 100];  % specification of contour values
    for k=1:length(tth) % k^th contour
        th = tth(k);
        Uzy = getFcontours(Uzy,Fzy,th,k); % adds contours to RGB image
        Uzx = getFcontours(Uzx,Fzx,th,k);
        Uyx = getFcontours(Uyx,Fyx,th,k);
    end
elseif conFLAG == 1
    
    tth = conVALUE; % specification of contour values
    Uzy = getFcontours(Uzy,Fzy,tth,5); % adds contours to RGB image
    Uzx = getFcontours(Uzx,Fzx,tth,5);
    Uyx = getFcontours(Uyx,Fyx,tth,5);
end

FF(1:Nz,wt+(1:Ny),:)                = uint8(Uzy);
FF(wt+Nz+(1:Nx),wt+(1:Nx),:)        = uint8(Uzx);
FF(wt+Nx+(1:Ny),2*wt+Nx+(1:Nx),:)	= uint8(Uyx);

%%% Tissue segmented
Tzx = squeeze(T(iy,:,:))'; % in z,x plane through source
Tzy = squeeze(T(:,ix,:))';
Tyx = squeeze(T(:,:,iz)); % in z,x plane through source
TT  = zeros(Nyy,Nxx);
TT(1:Nz,wt+(1:Ny))           	= (Tzy);
TT(wt+Nz +(1:Nx),wt+(1:Nx))     = (Tzx);
TT(wt+Nx+(1:Ny),2*wt+Nx+(1:Nx)) = (Tyx);
end

