
mcc -m OptogenSIM.m -a ./ -R '-startmsg, "Starting OptogenSIM V1.0: a 3D Monte Carlo simulation platform for light delivery design in optogentics"'

%%compile in windows 64
% mcc -m optogensim.m -a ./  -a mc_T.bin -a mc_VOL.bin -a gomcxyzOGS.exe -a mc_default.mci -a spectralLIB.mat ...
%     -R '-startmsg, "Starting OptogenSIM V1.0: a 3D Monte Carlo simulation platform for light delivery design in optogentics"'