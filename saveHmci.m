function simName = saveHmci(SIMname,myfolder,H,nm,tissue,originalname)
% function saveHmci(originalname,myname,myfolder,H,tissue)
% save mci file for mcxyzmay3
PRINTON = 0;

filename = sprintf('%s.mci',SIMname);
fprintf('saving %s ...',filename);

filename = fullfile(myfolder,filename);

PRINTON = 0;
if PRINTON
    disp(sprintf('-------- save %s --------',filename))
end

% tissue = makeTissueList_spmmouse1(nm); % also --> global tissue(1:Nt).s
tissue = makeTissueList_OGS(nm);
Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

if H(16) == inf % zfocus
    H(16) = 1e12;
end

fid = fopen(filename,'w');
% run parameters
fprintf(fid,'%s\n',SIMname);       % name (for name.mci, name_T.bin)
fprintf(fid,'%0.2f\n',H(1));    % time

fprintf(fid,'%d\n'   ,H(2));    % Nx
fprintf(fid,'%d\n'   ,H(3));    % Ny
fprintf(fid,'%d\n'   ,H(4));    % Nz
fprintf(fid,'%0.4f\n',H(5));    % dx
fprintf(fid,'%0.4f\n',H(6));    % dy 
fprintf(fid,'%0.4f\n',H(7));    % dz

% launch parameters
fprintf(fid,'%d\n'   ,H(8));    % mcflag
fprintf(fid,'%d\n'   ,H(9));  	% launchflag
fprintf(fid,'%d\n'   ,H(10));   % boundaryflag

fprintf(fid,'%0.4f\n',H(11));   % xs = photon launch position
fprintf(fid,'%0.4f\n',H(12));   % ys
fprintf(fid,'%0.4f\n',H(13));   % zs

fprintf(fid,'%0.4f\n',H(14)); % xfocus
fprintf(fid,'%0.4f\n',H(15)); % yfocus
fprintf(fid,'%0.4f\n',H(16)); % zfocus

% if launchflag==1, manually setting ux,uy,uz
fprintf(fid,'%0.4f\n',H(17)); % ux0
fprintf(fid,'%0.4f\n',H(18)); % uy0
fprintf(fid,'%0.4f\n',H(19)); % uz0

fprintf(fid,'%0.4f\n',H(20)); % radius
fprintf(fid,'%0.4f\n',H(21)); % waist
fprintf(fid,'%0.4f\n',H(22)); % wavelength
fprintf(fid,'%d\n',H(23));    % power [mW]

%%YL: add NA and thmax
fprintf(fid,'%0.4f\n',H(24)); % NA [-]
fprintf(fid,'%0.4f\n',H(25)); % fiber angle [radian]

fprintf(fid,'%d\n',H(26));    % Nt
fprintf(fid,'%d\n',H(27));    % YL: add timeFLAG

% tissue optical properties
for i=1:H(26)           % YL H(24)-> H(26) 
    fprintf(fid,'%0.4f\n',muav(i));
    fprintf(fid,'%0.4f\n',musv(i));
    fprintf(fid,'%0.4f\n',gv(i));
end
fclose(fid);

fprintf('saved.\n');

end