function [ix,iy,iz]=lookOGSoutput(F,VOL,T,bj,conFLAG,conVALUE,myfolder )
% function [ix,iy,iz]=lookbrainMay3c(F,VOL,T,dx,dy,dz,Nx,Ny,Nz,bj)
%   F   = fluence rate
%   VOL	= tissue image
%   T   = tissue cube (segmented)
%   bj  = initial choice of figure displayed
%       bj=1    tissue (gray)
%       bj=2    tissue (segmented)
%       bj=3    fluence contours on gray tissue
%       bj=4    fluence map
%       bj=5    return the check point to the source point 
%
% in ~/Dropbox/a_guestfolders/optogenSIM/optogenApr27 on steve's computer
% YL05262015: create a new figure to show the absolute fluence contour of
% desired light power density, passing light density(ld) FLAG and VALUE to control
%


s = '****************';
s = strvcat(s, '* lookbrain');
s = strvcat(s,'*****');
disp(s)

h1  = get(figure(500),'Position');
% pos = [h1(1)+h1(3)+20 h1(2) h1(3)*1.36 h1(4)-50];
pos = [h1(1)+h1(3)+20 h1(2) h1(3)*1.36 h1(4)-50];

clr = 'r';
mz = 12; sz = 12;

PRINTOUT = 0;
H = reportHmci('mc',PRINTOUT,myfolder);
Nx = H(2); Ny = H(3); Nz = H(4);
dx = H(5); dy = H(6); dz = H(7);
pwr = H(23);

F = F*pwr*10^-2; % YL: convert unit from 1/cm^2 to 1/mm^2

%% %%%%%%%%%%%%%%%
%% show brain
wt = 20; % width of text columns,rows around maps for x,y,z labels.

% photon source point
xs = H(11); 
ys = H(12);
zs = H(13);
%YL: view indices
ix = round(xs/dx + Nx/2+1/2);  
iy = round(ys/dy + Ny/2+1/2);  
iz = ceil(zs/dz+1);   


ixs = ix; iys = iy; izs = iz; % photon launch point
%YL,x,y,z
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1)*dz;

%YL,Tissue Image (RGB image)
Vzx = squeeze(VOL(iys,:,:))'; % in z,x plane through source
Vzy = squeeze(VOL(:,ixs,:))';
Vyx = squeeze(VOL(:,:,izs)); % in z,x plane through source

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
% Fzx = squeeze(F(iy,:,:))'; % in z,x plane through source
% Fzy = squeeze(F(:,ix,:))';
% Fyx = squeeze(F(:,:,iz)); % in z,x plane through source
%%YL
Fzx = squeeze(F(iys,:,:))'; % in z,x plane through source
Fzy = squeeze(F(:,ixs,:))';
Fyx = squeeze(F(:,:,izs)); % in z,x plane through source

FF  = zeros(Nyy,Nxx,3,'uint8');

%%YL: output
OUTallName = fullfile(myfolder,'OUTall.mat');
FzxCSVname = fullfile(myfolder,'Fzx_cp.csv');
FzxmapName = fullfile(myfolder,'Fzx_through_source.tif');

Fzx_zcp = Fzx(:,ixs);
save(OUTallName,'Fzx','Fyx','Fzy','Fzx_zcp','ixs','iys','izs','xs','ys','zs','x','y','z');
csvwrite(FzxCSVname,[z' Fzx_zcp]);
fig503 = figure(503); set(fig503,'Position',[100 100 size(Fzx)],'Visible','off');
imagesc(x,z,Fzx),colormap('jet'),colorbar,xlabel('X, [cm]'),ylabel('Z, [cm]'),
xlim([xs-0.10 xs+0.10]); ylim([max(0,zs-0.10) zs+0.15])
print(fig503,'-dtiff','-r300',FzxmapName)

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

%% show brain
% buttons
bh = 10; bv = bh;
bw = 2*bh;
% blabel(1).s = 'tissue (gray)';
% blabel(2).s = 'tissue (color)';
% blabel(3).s = 'segmented';
% blabel(4).s = 'fluence contours';
% blabel(5).s = 'fluence';
blabel(1).s = 'tissue (gray)';
blabel(2).s = 'segmented';
blabel(3).s = 'fluence contours';
blabel(4).s = 'fluence';
blabel(5).s = '@source';

Nb = length(blabel);
bx = (zeros(1,Nb)+1)*bw;
by = Nyy - [Nb:-1:1]*(bw+8);

ss='to RETURN to control --> click outside figure,hit Return.';

flag = 1;
flagNOCHANGE = 0;

fig = figure(501);clf
set(fig,'Position',pos)
set(fig,'name','OptogenSIM output','numbertitle','off')


%%
while flag % CYCLE within this loop until click outside figure.
    %try % check for errors. If so, disp error messange, then retry.
        
        sss = ss; % restore original title
        
        % FIGURE
        figure(fig); hold off
        switch bj
            case 1
                imagesc(VVOL)
           
            case 2
                imagesc(TT)
                hold on
                colormap(makec2f)
                rectangle('position',[320 700 20 20],'facecolor','y')
                text(350,710,'white matter','color','y','fontsize',12)
                rectangle('position',[320 750 20 20],'facecolor',[1 .7 0])
                text(350,760,'gray matter','color',[1 .7 0],'fontsize',12)
            case 3
                imagesc(FF)
                hold on
                xx = Nz*0.99; yy = Nyy-70; dyy = 20;
                %YL: unit conversion from cm to mm
                if conFLAG == 0   %YL: build-in four contours
                    text(xx,yy-4*dyy,'0.1 mW/mm^{2}','color','w','fontsize',10)
                    text(xx,yy-3*dyy,'1   mW/mm^{2}','color','g','fontsize',10)
                    text(xx,yy-2*dyy,'10  mW/mm^{2}','color','y','fontsize',10)
                    text(xx,yy-1*dyy,'100 mW/mm^{2}','color','r','fontsize',10)
                    text(xx-30,Nyy-30,sprintf('Power = %4.3f mW',pwr),'color','w','fontsize',14)
                elseif conFLAG == 1  %YL: specified contour
                    text(xx,yy-2.5*dyy,sprintf('%4s mW/mm^{2}',num2str(conVALUE)),'color','m','fontsize',12)
                    text(xx-30,Nyy-30,sprintf('Power = %4.3f mW',pwr),'color','w','fontsize',14)
                end
            case 4
                %YL
                flow = 0.1; fhigh = 100;
                ind = find(FFF >= flow);
                xx = squeeze(VVOL(:,:,1));
%                 xx(ind) = 0;
                FFF2 = FFF;
                FFF2(ind) = FFF2(ind)+10*double(max(xx(:)));
                FFF2(~ind) = NaN;
                CC = imfuse(xx,FFF2,'falsecolor');
                imagesc(CC);
%                 imshowpair(xx,FFF,'falsecolor')
      
        end
        hold on    
        
        % plot crosses on maps
        clr = 'r'; lw = 2;
        plot([1 1]*(wt+iys),[0 izs],[clr '-'],'markersize',mz,'linewidth',lw)
        plot([1 1]*(wt+ixs),wt+Nz+[0 izs],[clr '-'],'markersize',mz,'linewidth',lw)
        clr = 'c';
        plot((wt+iy),iz,[clr '+'],'markersize',mz,'linewidth',lw)
        plot((wt+ix),wt+Nz+iz,[clr '+'],'markersize',mz,'linewidth',lw)
        plot(2*wt+Nx+ix,wt+Nz+iy,[clr '+'],'markersize',mz,'linewidth',lw)

        % xyz labels
        sz = 12;
        text(5,Nz/2,'z','color','w','fontsize',sz)
        text(Ny/2,Nz+wt-10,'y','color','w','fontsize',sz)

        text(5,wt+Nz+Nz/2,sprintf('z'),'color','w','fontsize',sz)
        text(wt+Nx/2,wt+2*Nz+wt+10-20,'x','color','w','fontsize',sz)

        text(wt+Nx+7,Nyy-Ny/2-Ny/5,'y','color','w','fontsize',sz)
        text(wt+Nx+Nx/2,2*wt+Nz+Ny-10,'x','color','w','fontsize',sz)
        
        plot(wt+[1 1 Ny],[1 Nz Nz],'w-')
        plot(wt+[1 1 Nx],wt+Nz+[1 Nz Nz],'w-')
        plot(2*wt+Nz+[1 1 Nx],wt+Nz+[1 Ny Ny],'w-')
        
        % show xs,ys,zs values
        xx = 5;
        yy = 2*(Nz+wt)+20;
     	text(xx,yy,sprintf('x_s = %0.3f cm,\tx = %0.3f cm',xs,x(ix)),'color','w','fontsize',12)
     	text(xx,yy+30,sprintf('y_s = %0.3f cm,\ty = %0.3f cm',ys,y(iy)),'color','w','fontsize',12)
     	text(xx,yy+60,sprintf('z_s = %0.3f cm,\tz = %0.3f cm',zs,z(iz)),'color','w','fontsize',12)
        xx = 20;
%      	text(xx,yy+90,sprintf('%s = %0.2e mW/cm^{2}','\phi',F(iy,ix,iz)),'color','w','fontsize',14)
     	text(xx,yy+90,sprintf('%s = %0.2e mW/mm^{2}','\phi',F(iy,ix,iz)),'color','w','fontsize',14) % YL: unit conversion
            

        % buttons
        for j=1:Nb
            rectangle('position',[bx(j)-bh,by(j)-bv 2*bh 2*bv],'edgecolor','w','facecolor','k')
            text(bx(j)+2*bh,by(j),blabel(j).s,'color','w','fontsize',12)
            if j==bj
                plot(bx(bj),by(bj),'wx','markersize',12,'linewidth',2)
            end
            %YL
           if j == 5 & ix == ixs & iy == iys & iz == izs
               plot(bx(j),by(j),'wx','markersize',12,'linewidth',2)
           end
            
            
        end
     
        if flagNOCHANGE
            sss = 'no change. try again.';
            flagNOCHANGE = 0;
        end
        
        title(sss,'fontsize',12,'fontname','courier')
        axis equal image off

        % CLICK
        flag = 0;
        maxerror = 5; errnum = 0;
        while flag==0
            try
                figure(fig)
                uu = round(ginput); u = uu(end,:); % only use last click
                flag = 1;
            catch ME
                disp('oops. try again.')
                errnum = errnum + 1;
                if errnum == maxerror
                    close(figure(501))
                    break
                end
                
            end
        end
        
        h  = u(1); % integer values
        v  = u(2);
        hv = uu(1); % double values (not used yet)
        vv = uu(2);
        if h<=0 | h>Nxx | v<=0 | v>Nyy % clicked outside figure = quit
            flag=0;
            title('')
            disp('done')
            ixyzflag =0; % YL: flag for the change of ix, iy, iz
        elseif v >= by(1)-bv & v<= by(1)+bv & h>=bx(1)-bh & h<=bx(1)+bh % buttons
            bj = 1;
            ixyzflag =0; % YL: flag for the change of ix, iy, iz
        elseif v >= by(2)-bv & v<= by(2)+bv & h>=bx(2)-bh & h<=bx(2)+bh
            bj = 2;
            ixyzflag =0; % YL: flag for the change of ix, iy, iz
        elseif v >= by(3)-bv & v<= by(3)+bv & h>=bx(3)-bh & h<=bx(3)+bh
            bj = 3;            
            ixyzflag =0; % YL: flag for the change of ix, iy, iz
        elseif v >= by(4)-bv & v<= by(4)+bv & h>=bx(4)-bh & h<=bx(4)+bh
            bj = 4;
            ixyzflag =0; % YL: flag for the change of ix, iy, iz
        elseif v >= by(5)-bv & v<= by(5)+bv & h>=bx(5)-bh & h<=bx(5)+bh
            bj = bj;   % YL: do not change bj,just chang ix,iy,iz to ixs, iys,izs 
            ixyzflag =1; % YL: flag for the change of ix, iy, iz
            ix = ixs; iy = iys; iz = izs;
        elseif v<=Nz & h>wt & h<wt+Ny % Vzy, using integer pointers. Better to interpolate.
                iz = v;
                iy = h-wt;
                %zv = z(v);
                %yv = y(h);
                disp('Last click is on ZY section')
                ixyzflag =1; % YL: flag for the change of ix, iy, iz
        elseif h>wt & h<Nx+wt & v<2*Nz+wt & v> Nz+wt % Vzx
                iz = v-Nz-wt;
                ix = h-wt;            
                %zv = z(v-Nz);
                %xv = x(h);
               disp('Last click is on ZX section')
                ixyzflag =1; % YL: flag for the change of ix, iy, iz
        elseif h>2*wt+Nx & v>Nz+wt & v<Nz+Ny+wt% Vyx
                iy = v-Nz-wt;
                ix = h-Nx-2*wt;            
                %yv = y(v-Nz);
                %xv = x(h-Nx);
                disp('Last click is on YX section')
                ixyzflag =1; % YL: flag for the change of ix, iy, iz
        else 
                disp('no change')
                flagNOCHANGE = 1;
                ixyzflag =0; % YL: flag for the change of ix, iy, iz
        end
        
        %YL: update the atlas sections at the positions other than the lauching position  
       if ixyzflag ~=0
        if ix ~= ixs | iy ~= iys | iz ~= izs
            
%             disp(sprintf('Loading the sections determined by ix = %6.4f, iy = %6.4f, iz = %6.4f',ix,iy,iz))
            [VVOL TT FF FFF] = getVFTxyz(F,VOL,T,conFLAG,conVALUE,ix,iy,iz,Nx,Ny,Nz,dx,dy,dz,wt);
            
        elseif ix == ixs & iy == iys & iz == izs
            disp(sprintf('Returning to the photon launching position: ix = %6.4f, iy = %6.4f, iz = %6.4f',ix,iy,iz))
            [VVOL TT FF FFF] = getVFTxyz(F,VOL,T,conFLAG,conVALUE,ix,iy,iz,Nx,Ny,Nz,dx,dy,dz,wt);
      
        end
        
        
       end
        
        
    %%
%     catch
%         disp('error')
%         for i=1:3
%             sss = 'oops. try again.';
%         end
%         keyboard
%     end
end % While flag

% %% print:figure501, figure1, and figure 500
% sz0 = get(0,'screensize');
% sw0 = sz0(3);
% sh0 = sz0(4);
% 
% RES = 300;
% figure(500); pixw = sw0*0.25;pixh = sh0*0.85; set(gcf,'PaperUnits','inches','PaperPosition',[0 0 pixw/RES pixh/RES]);filename = 'GUI.jpg'; print('-djpeg',['-r',num2str(RES)],filename)


