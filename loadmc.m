function [T VOL H SIMname] = loadmc(SIMname,myfolder)
% loads mc.mci, mc_T.bin, mc_VOL.bin --> T,VOL,H
% presumes they exist.


% load mc.mci
PRINTON = 0;
[H SIMname] = reportHmci(SIMname,PRINTON,myfolder);
Nx = H(2);
Ny = H(3);
Nz = H(4);

% Load tissue structure in voxels, T(y,x,z) 
% filename = sprintf('mcLibrary/%s_T.bin',SIMname);
% filename = fullfile(pwd,'mcLibrary',sprintf('%s_T.bin','mc'));  % YL: tentatively only for mc_T.bin 
filename = fullfile(myfolder,sprintf('%s_T.bin','mc'));  % YL: tentatively only for mc_T.bin 

disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)
clear Data

% flagV = (exist('mcLibrary/mc_VOL.bin')>0);
flagV = (exist('mc_VOL.bin')>0);

if flagV
    % Load tissue structure in voxels, VOL(y,x,z) 
%     filename = sprintf('mcLibrary/%s_VOL.bin',SIMname);
%     filename = fullfile(pwd,'mcLibrary',sprintf('%s_VOL.bin','mc'));% YL: tentatively only for mc_T.bin 
    filename = fullfile(myfolder, sprintf('%s_VOL.bin','mc'));% YL: tentatively only for mc_T.bin 

    disp(['loading ' filename])
    tic
        fid = fopen(filename, 'rb');
        [Data count] = fread(fid, Ny*Nx*Nz, 'uint16');
        fclose(fid);
    toc
    VOL = reshape(Data,Ny,Nx,Nz); % VOL(y,x,z)
    clear Data
else
    fprintf('No VOL. Let VOL = T.\n')
    VOL = double(T);
end
