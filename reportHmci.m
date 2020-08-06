function [H SIMname2] = reportHmci(SIMname,PRINTON,myfolder)
% function [H originalname] = reportHmci(myname,PRINTON)
%   Lists the values of the input file  myname.mci.
%   used by mcxyz_may3.c.
if nargin==1
    PRINTON= 1;
end
% filename = sprintf('mcLibrary/%s.mci',SIMname);
% filename = fullfile(pwd,'mcLibrary',sprintf('%s.mci',SIMname)); % YL
filename = fullfile(myfolder,sprintf('%s.mci',SIMname)); % YL

fid = fopen(filename,'r');
SIMname2 = fgetl(fid);
H = fscanf(fid,'%f');
fclose(fid);

s(1).s = 'time_photons';
s(2).s = 'Nx';
s(3).s = 'Ny';
s(4).s = 'Nz';
s(5).s = 'dx';
s(6).s = 'dy';
s(7).s = 'dz';
s(8).s = 'mcflag';
s(9).s = 'launchflag';
s(10).s = 'boundaryflag';
s(11).s = 'xs';
s(12).s = 'ys';
s(13).s = 'zs';
s(14).s = 'xfocus';
s(15).s = 'yfocus';
s(16).s = 'zfocus';
s(17).s = 'ux0';
s(18).s = 'uy0';
s(19).s = 'uz0';
s(20).s = 'radius';
s(21).s = 'waist';
s(22).s = 'wavelength';
s(23).s = 'power [mW]';
s(24).s = 'fiber NA [-]'; %YL
s(25).s = 'fiber half angle[radian]'; %YL
s(26).s = '# of tissues';
s(27).s = '[1: specify simulation time; 2: speicfy photons]';  %YL

if PRINTON
    fprintf('SIM name = %s\n',SIMname2);
    for i=1:27    %YL: 23 --> 27
        disp(sprintf('%d\t%10s = %0.4f',i,s(i).s,H(i)))
    end
    
    for j=1:H(26)  % YL: H(24) -> H(26)
        i=i+1;
        disp(sprintf('---'))
        disp(sprintf('%d\tmua = %0.4f',1i,H(i)))
        i=i+1;
        disp(sprintf('%d\tmus = %0.4f',i,H(i)))
        i=i+1;
        disp(sprintf('%d\tg   = %0.4f',i,H(i)))
    end
end

