function tissue = makeTissueList_OGS(nm)
% function tissue = makeTissueList_Azimipour(nm)
% 4-29-2015: use Azimipour optical properties for Gray and White matter.
% 12-16-2014: use the best estimates for gray matter and white matter by Steven L. Jacques. 
% 12-26-2014: correct the B, S of gray and white matter
%
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%   Uses 
%       SpectralLIB.mat
%
% Steven L. Jacques, March 3, 2013.

PRINTON = 0;

%% Load spectral library
load SpectralLIB
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muamel,nm);
LOADED = 1;


%% Load interpolated brain tissue optical properties for wavelengh of 473, 594, and 635nm
%load('Interpolated_brain_optical_properties.mat','BOP');
% BOP = 
%     wavelengh: [473 594 635]
%        grayus: [112.2976 94.7725 89.2915]
%        grayua: [0.5642 0.2206 0.1988]
%         grayG: [0.8800 0.8843 0.8910]
%       whiteus: [424.1268 415.9300 408.0428]
%       whiteua: [1.1990 0.8544 0.7883]
%        whiteG: [0.7936 0.8321 0.8412]

%% Create tissueList

j=1;
tissue(j).name = 'air';
tissue(j).mua = 0.1;
tissue(j).mus = 10.0;
tissue(j).g   = 1.0;

j=2;
tissue(j).name = 'water';
tissue(j).mua = 0.1;
tissue(j).mus = 10.0;
tissue(j).g   = 1.0;

j=3;
tissue(j).name = 'blood';
B       = 1.00;
S       = 0.75;
W       = 0.95;
M       = 0;
musp500 = 10;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j = 4;
tissue(j).name = 'dermis';
B = 0.002; 
S = 0.67;
W = 0.65;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=5;
tissue(j).name = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=6;
tissue(j).name = 'skull';
B = 0.0005;
S = 0.75;
W = 0.35;
M = 0;
musp500 = 30;
fray    = 0.0;
bmie    = 1.0;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=7;
tissue(j).name = 'gray matter';
% B = 0.037;
% S = 0.65;%;
B = 0.028;  %YL 12-26: B = 2.763 +/- 0.669
S = 0.62;   %YL 12-26: S = 62.100 +/- 5.340
W = 0.65;%
M = 0;
musp500 = 23.7;  % 22.2±2.562  [cm^-1]
fray    = 0;
bmie    = 1.15;    % 1.40±0.261
X       = [B*S B*(1-S) W M]';
musp    = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg      = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=8;
tissue(j).name = 'white matter';
B = 0.028;  %YL 12-26: B = 2.763 +/- 0.669
S = 0.62;   %YL 12-26: S = 62.100 +/- 5.340
W = 0.65;%
M = 0;
musp500 = 50.5;  % 
fray    = 0.0;
bmie    = 0.7;
X = [B*S B*(1-S) W M]';   
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
gg = 0.90;
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

%% CSF
j=9;
tissue(j).name = 'CSF';
B       = 0;
S       = 0;
W       = 1;
M       = 0;
musp500 = 2.4;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
X = [B*S B*(1-S) W M]';
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

Nt = length(tissue);

for j=1:Nt, tissue(j).nm = nm; end % wavelength

if PRINTON
disp(sprintf('---- tissueList ------ \t%10s   \t%s  \t%s  \t%s','mua','mus','g','musp'))
for i=1:Nt
    musp = tissue(i).mus*(1-tissue(i).g);
    disp(sprintf('%d\t%15s\t%10.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,musp))
end
disp(' ')
end

