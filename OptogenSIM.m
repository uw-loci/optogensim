function OptogenSIM

% OptogenSIM.m
% This Matlab GUI is associated with a 3D  Monte Carlo simulation platform
% in optogentics application, used for simulating light transport in rodent
% brain tissue atlases
% This platform is a collaborating work being conducted among three groups
% from Oregon Health and Science University, University of Wisconsin at
% Milwaukee, and University of Wisconsin at Madison. Please check the
% website at http://loci.wisc.edu/software/optogensim for details.
%Main developers:
%Yuming Liu, Steven.L. Jacques
%Laboratory for Optical and Computational Instrumentation
%University of Wisconsin-Madison
%Since Sept, 2014
%August 2015: verison1.0 is based on "optogensimMay15slj"

% To deploy this code, use:
% mcc -m optogensim.m -a ./ -R '-startmsg, "Starting OptogenSIM V1.0: a 3D Monte Carlo simulation platform for light delivery design in optogentics"'


clear all; home;

if (~isdeployed)
    addpath('./SliceBrowser');
    addpath('.');
end
%%
COMP = ismac;  % 0 = Windows, 1 = Mac
%%
%% locate the mc_VOL.bin, mc_T.bin,mc_default.mci directory

% if (~isdeployed)
%     
%     [ctfpath VOLname] = fileparts(which('mc_VOL.bin'));
% elseif (isdeployed)
% %     ctfpath = ctfroot;
%     ctfpath = pwd;
% 
% end

[ctfpath VOLname] = fileparts(which('mc_VOL.bin'));


% if COMP     %  mac
%     godir = './';
%     librdir = fullfile(ctfpath,'OGSworking');
%     
% else         % pc
%     godir = './';
%     librdir = fullfile(ctfpath,'OGSworking');
% end

librdir = fullfile(pwd,'OGSworking');
if ~exist(librdir,'dir')
  mkdir(librdir)
end


if ~exist(fullfile(librdir,'mc.mci'),'file')
    
    copyfile(fullfile(ctfpath,'mc_default.mci'), fullfile(librdir,'mc.mci'))
 
end

if ~exist(fullfile(librdir,'mc_VOL.bin'),'file')
    
    copyfile(fullfile(ctfpath,'mc_VOL.bin'), fullfile(librdir,'mc_VOL.bin'))
 
end

if ~exist(fullfile(librdir,'mc_T.bin'),'file')
    
    copyfile(fullfile(ctfpath,'mc_T.bin'), fullfile(librdir,'mc_T.bin'))
 
end

if COMP
    if ~exist(fullfile(librdir,'gomcxyzOGS'),'file')
        copyfile(fullfile(ctfpath,'gomcxyzOGS'), fullfile(librdir,'gomcxyzOGS'))
    end
    
else
    if ~exist(fullfile(librdir,'gomcxyzOGS.exe'),'file')
        copyfile(fullfile(ctfpath,'gomcxyzOGS.exe'), fullfile(librdir,'gomcxyzOGS.exe'))
    end
    
end


% default parameters
% refn        = 1.36;     % tissue refractive index
timeFLAG = 1;      % 1: specify simulation time or photons; 2: specify the number of photons in simulation 
simFLAG     = 1  ; % 1: default or previous simulation parameters; 0: updated simulations
flagOUTGO   = 1; % 1 = on, 0 = off
conVALUE = 5;  % 5mw/mm-2, YL: default value for the specified contour, 
conFLAG = 0; % YL: 0: use built-in four contours; 1: use the specified contour
VOL2 = [];  % YL
T   = [];
VOL = [];
Nx  = [];
Ny  = [];
Nz = [];
dx = [];
dy = [];
dz = [];
xs = [];
ys = [];
zs = [];
mcflag = [];  
mcflagtxt = {'De-focused Gaussian(1/e^2 radius)','De-focused Gaussian(1/e radius)','Defocused flat beam',...
    'Isotropic point source','Focused Gaussian(1/e^2 radius)','Focused Gaussian(1/e radius)','Focused flat beam'};
currentmcflag = [];
originalname ='';

SIMname = 'mc';
[T VOL H SIMname] = loadmc(SIMname,librdir); % Loads last mc simulation (mc.mci, mc_T.bin, mc_VOL.bin)
if H(16) ==  1.0000e+12
    H(16) = inf;
end
time_photons = H(1); % time or photons  of simulation
Nx      = H(2);
Ny      = H(3);
Nz      = H(4);
dx      = H(5);
dy      = H(6);
dz      = H(7);
mcflag  = H(8);   %0: focused flat; 1:focused Gaussian 1/e; 2: isotropic point; 3:defocused flat ...
                  %4: defocused Gaussian 1/e; 5: focused Gaussian 1/e2; 6: de-focused Gaussian 1/e2;
launchflag = H(9);
boundaryflag = H(10);
xs      = H(11);
ys      = H(12);
zs      = H(13);
xfocus  = H(14);
yfocus  = H(15);
zfocus  = H(16);
ux0     = H(17);
uy0     = H(18);
uz0     = H(19);
radius  = H(20);
waist   = H(21);
nm = H(22);    % light wavelength
pwr     = H(23);       % power of the light source
NA =     H(24);        % numeric aperture of the fiber
thmax =  H(25);
Nt      = H(26);
timeFLAG = H(27);      % 1: specify simulation time or photons; 2: specify the number of photons in simulation 

tissue      = struct('name',[],'mua',[],'mus',[],'g',[]);
tissue      = makeTissueList_OGS(nm);

%%%%%%%%%%%
% control panel
%%%
sz0 = get(0,'screensize');
sw0 = sz0(3);
sh0 = sz0(4);
whr = sw0/sh0; % screen width / height ratio

% scale = [sw0 sh0]*1.5;
% pg = [0.2 0.25 scale];   % position of control gui

% scale = [.4 1]*0.9;
% pg = [0.15 0.05 scale];   % position of window
%YL
pg = [0.01 0.10 0.25 0.85];   % position of control gui

xB = 0; xB2 = 0.5; yB=0.4; yB2 = 0.35; dyB = .05; % position of button
wB = .5; hB = .1; hB2 = 0.05; % size of button

% CONTROL PANEL
guiCtrl = figure(500);clf
set(guiCtrl,'Resize','on','Units','pixels','Position',[sw0*pg(1) sh0*pg(2)  sw0*pg(3) sh0*pg(4)],'Visible','off',...
    'MenuBar','none','name','OptogenSIM V1.0','NumberTitle','off','UserData',0);

%%%%%%%%%%%%
% Edit boxes for simulation settings
%%%%
xB = 0.03; xB2 = 0.35; yB=1; dyB = .035; % position of labels & textboxes
wB = .25; hB = .03;  % size of button
fz = .5;

% Edit box for SIMname
j=1;
listtextbox_SIM = uicontrol('Parent',guiCtrl,'Style','text','String','Simulation name',...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterSIM = uicontrol('Parent',guiCtrl,'Style','edit','String',SIMname,'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_SIM});


% Edit box for power
j=j+1;
listtextbox_pwr = uicontrol('Parent',guiCtrl,'Style','text','String','Power [mw]',...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterpwr = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(pwr),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_pwr});

% Edit box for wavelength
j=j+1;
listtextbox_nm = uicontrol('Parent',guiCtrl,'Style','text','String',['Wavelength [nm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterNM = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(nm),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_nm});

% Edit box for time or photons
j=j+1;
min_numChoice = uicontrol('Parent',guiCtrl,'Style','popupmenu','String',{'Time [min]','Photons [#]'},...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB],'Value',timeFLAG,'Callback',{@timeorphoton_callback});

entertime_photon = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(time_photons),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_min});

% Edit box for radius
j=j+1;
listtextbox_radius = uicontrol('Parent',guiCtrl,'Style','text','String',['Fiber radius [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterradius = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(radius),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_radius});
%edit box for fiber NA
% Edit box for radius
j=j+1;
listtextbox_NA = uicontrol('Parent',guiCtrl,'Style','text','String',['Fiber NA [-]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterNA = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(NA),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_NA});


% Edit box for xs
j=j+1;
listtextbox_xs = uicontrol('Parent',guiCtrl,'Style','text','String',['Source X [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterxs = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(xs),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_xs});
% Edit box for ys
j=j+1;
listtextbox_ys = uicontrol('Parent',guiCtrl,'Style','text','String',['Source Y [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterys = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(ys),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_ys});
% Edit box for zs
j=j+1;
listtextbox_zs = uicontrol('Parent',guiCtrl,'Style','text','String',['Source Z [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterzs = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(zs),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_zs});

% Edit box for mcflag
j=j+1;
listtextbox_mcflag = uicontrol('Parent',guiCtrl,'Style','text','String',['Beam type'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB*0.5 yB-dyB*j*0.975 wB*0.75 hB]);

setmcflag(mcflag)  % update currentmcflag
beamchoice = uicontrol('Parent',guiCtrl,'Style','popupmenu','String',mcflagtxt,'BackgroundColor','w',...
    'Min',0,'Max',1,'Value',currentmcflag,'UserData',[],'Units','normalized','Position',[xB yB-dyB*(j+1.25) wB*2.25 hB*2],'Callback',{@beamchoice_callback});

% beamchoice = uicontrol('Parent',guiCtrl,'Style','edit','String',currentmcflag,'BackgroundColor','w',...
%     'Min',0,'Max',1,'UserData',[],'enable','off','Units','normalized','Position',[xB2 yB-dyB*(j+1.25) wB hB],'Callback',{@beamchoice_callback});

% Edit box for xfocus
j=j+1.5;
listtextbox_xfocus = uicontrol('Parent',guiCtrl,'Style','text','String',['Focus X [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterxfocus = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(xfocus),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_xfocus});
% Edit box for yfocus
j=j+1;
listtextbox_yfocus = uicontrol('Parent',guiCtrl,'Style','text','String',['Focus Y [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enteryfocus = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(yfocus),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_yfocus});
% Edit box for zfocus
j=j+1;
listtextbox_zfocus = uicontrol('Parent',guiCtrl,'Style','text','String',['Focus Z [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterzfocus = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(zfocus),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_zfocus});


% Edit box for waist
j=j+1;
listtextbox_waist = uicontrol('Parent',guiCtrl,'Style','text','String',['Waist [cm]'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterwaist = uicontrol('Parent',guiCtrl,'Style','edit','String',num2str(waist),'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_waist});

j = j+1.0;
% TEXT message
remindLabel = uicontrol('Parent',guiCtrl,'Style','text','String','Hit return after changes',...
    'FontUnits','normalized','FontSize',.50,'Units','normalized','Position',[xB2*0.8 yB-dyB*j wB*1.4 hB]);

% Edit box for desired light density
j=j+0.5;
textbox_ld1 = uicontrol('Parent',guiCtrl,'Style','text','String','Flunece contour(s)',...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB-0.075 yB-dyB*(j+0.75) wB+0.15 hB]);
textbox_ld2 = uicontrol('Parent',guiCtrl,'Style','text','String','[mw/mm^-2]',...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB-0.075 yB-dyB*(j+1.25) wB+0.15 hB]);

% enterlightden = uicontrol('Parent',guiCtrl,'Style','list','String',num2str(3),'BackgroundColor','w',...
%     'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB hB],'Callback',{@textbox_ligthDEN});
% ldcontours = [{'0.1'};{'1'};{'10'};{'100'}];
defaultcontours = sprintf('%4s %4s %4s %4s','0.1', '1','10','100');
lightdenmenu = uicontrol('Parent',guiCtrl,'Style','popupmenu','String',[{'Built-in 4 Contours'};{'Specify a contour'}],'BackgroundColor','w',...
    'Min',0,'Max',1,'Value',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*(j+1.25) wB*1.5 hB*2],'Callback',{@setcontour_callback});
enterlightden = uicontrol('Parent',guiCtrl,'Style','edit','String',defaultcontours,'BackgroundColor','w',...
    'Min',0,'Max',1,'UserData',[],'enable','off','Units','normalized','Position',[xB2 yB-dyB*(j+1.25) wB*1.5 hB],'Callback',{@entercontour_callback});

%%YL: control buttons

xB3 = 0.655; dyB3 = 0.050; wB3 = 0.30; hB3 = 0.04; fz3 = 0.30;
loadsimulation = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Load SIM',...
    'FontUnits','normalized','FontSize',fz3,'Fontweight','bold','Units','normalized',...
    'Position',[xB3 yB-dyB3*1.5 wB3 hB3*1.5],'callback','ClickedCallback','Callback',{@loadSIM});
load3datlas = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Load Atlas',...
    'FontUnits','normalized','FontSize',fz3,'Fontweight','bold','Units','normalized',...
    'Position',[xB3 yB-dyB3*1.5*2 wB3 hB3*1.5],'callback','ClickedCallback','Callback',{@loadATLAS});
loadoplib = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Brain OPLIB',...
    'FontUnits','normalized','FontSize',fz3,'Fontweight','bold','Units','normalized',...
    'Position',[xB3 yB-dyB3*1.5*3 wB3 hB3*1.5],'callback','ClickedCallback','Callback',{@loadOPLIB});

% % SAVE BUTTON to view simulation parameters
callsavemci = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Save SIM','BackgroundColor','g',...
    'FontUnits','normalized','FontSize',fz3,'fontweight','bold','Units','normalized','Position',[xB3 yB-dyB3*1.5*4 wB3 hB3*1.5],...
    'Callback',{@savemci});


% % CHECKMCI BUTTON to view mc.mci
callcheckmci = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Check SIM',...
    'FontUnits','normalized','FontSize',fz3,'fontweight','bold','Units','normalized','Position',[xB3 yB-dyB3*1.5*5 wB3 hB3*1.5],...
    'Callback',{@checkmci});

% GO BUTTON to run mcxyz
GOcolor = 'g';
callGO = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','G O','BackgroundColor','g',...
    'FontUnits','normalized','FontSize',fz3,'fontweight','bold','Units','normalized','Position',[xB3 yB-dyB3*1.5*6 wB3 hB3*1.5],...
    'Callback',{@GO});

callOUTPUT = uicontrol('Parent',guiCtrl,'Style','pushbutton','String','Check LightDen.','BackgroundColor','g',...
    'FontUnits','normalized','FontSize',fz3,'fontweight','bold','Units','normalized','Position',[xB3 yB-dyB3*1.5*7 wB3 hB3*1.5],...
    'Callback',{@OUTPUT});


%%%% MESSAGES
% % TEXT message
% remindLabel = uicontrol('Parent',guiCtrl,'Style','text','String','hit return after changes',...
%     'FontUnits','normalized','FontSize',.30,'Units','normalized','Position',[0.65 .9 .25 .06]);

%YL
% defaultBackground = get(0,'defaultUicontrolBackgroundColor');
% set(guiCtrl,'Color',defaultBackground);
% set(guiCtrl,'Visible','on');

% % MESSAGE box
% note1 = 'ready';
% infoLabel = uicontrol('Parent',guiCtrl,'Style','text','String',note1,'FontUnits',...
%     'normalized','FontSize',.068,'Units','normalized','Position',[0 .0 1.0 .40],'BackgroundColor',[.7 1 1]);

%YL: % BUTTON for DEBUG
% note2 = 'turn debug on';
% calldebugON = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',note2,...
%     'FontUnits','normalized','FontSize',.30,'Fontweight','bold','Units','normalized',...
%     'Position',[.80 .05 .20 .05],'callback','ClickedCallback','Callback',{@debugON});



% set(infoLabel,'String','Reading files.');

%%%% Load mc_H.mci and mc_T.bin
%  	 Assumes mc_T.bin and mc_H.mci exist. If not, quits.
%       exist() yields 2 if file exists. 0 if not.
flagH = (exist(fullfile(librdir,'mc.mci'))>0);
flagT = (exist(fullfile(librdir,'mc_T.bin'))>0);
flagV = (exist(fullfile(librdir,'mc_VOL.bin'))>0);
flagF = (exist(fullfile(librdir,'mc_F.bin'))>0);
if ~flagH, disp('missing mc.mci');end
if ~flagT, disp('missing mc_T.bin. Load a simulation.');end
if ~flagV, disp('missing mc_VOL.bin. Uses mc_T.bin instead.');end
if ~flagF, 
    disp('missing mc_F.bin. Be sure to run GO before using OUTPUT.'); 
    set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
    set(callOUTPUT,'Enable','off')
end


%%YL: load the atlas to check/update the launching position
VOL2 = permute(VOL,[2 1 3]);
xe = round(xs/dx+1/2+Nx/2);
ye= round(ys/dy+1/2 + Ny/2);
ze = round(zs/dz+1/2);
[sbxyz sbfig] = SliceBrowser(VOL2,xe,ye,ze); % Nx, Ny, Nz
sbfigh = sbfig; %figure(sbfig);   % handle to the sbfig
% disp(sprintf('%s',sbxyz))

%YL:add  a push putton on the atals to check/update the launching position
enterxyz = uicontrol('Parent',sbfig,'Style','pushbutton','String','Use current XYZ as launching position','Visible','on',...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized',...
    'Position',[0.468 .066 .398 .054],'callback','ClickedCallback','Callback', {@enterXYZ_Callback});


if ~isempty(find(currentmcflag == [1 2 3 4]))
     set([enterxfocus enteryfocus enterzfocus enterwaist], 'String',num2str([]),'Enable','off')
     H(14) =NaN; H(15) = NaN; H(16) = NaN; H(21) = NaN; 
else
    set([enterxfocus enteryfocus enterzfocus enterwaist],'Enable','on')
end

% MESSAGE box
note1 = sprintf('Ready');
disp(sprintf('Current working directory is %s',librdir))
infoLabel = uicontrol('Parent',guiCtrl,'Style','text','String',note1,'FontUnits',...
    'normalized','fontname','courier','fontsize',0.08,'Units','normalized','Position',[0 .0 1.0 .40],'BackgroundColor',[.7 1 1]);


% BUTTON to reset gui
note3 = 'Reset';
imgReset = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',note3,...
    'FontUnits','normalized','FontSize',.50,'Fontweight','bold','Units','normalized',...
    'Position',[.80 .0 .20 .05],'callback','ClickedCallback','Callback',{@resetOGS});

defaultBackground = get(0,'defaultUicontrolBackgroundColor');
set(guiCtrl,'Color',defaultBackground);
set(guiCtrl,'Visible','on');

%%%%%%%%%%%%%%read
%% CALLBACKS
%--------------------------------------------------------------------------
    function loadOPLIB(loadoplib,eventdata)
        set(infoLabel,'String','The current optical properties library of brain tisse is included in the file: makeTissueList_OGS.m');
    end

%--------------------------------------------------------------------------
%% set beam choice
    function beamchoice_callback(beammenu,eventdata)
         
        currentmcflag = get(beamchoice,'Value');
        if  currentmcflag == 1
            mcflag = 6;
        elseif currentmcflag == 2
            mcflag = 4;
        elseif currentmcflag == 3
            mcflag = 3;
        elseif currentmcflag == 4
            mcflag = 2;
        elseif currentmcflag == 5
            mcflag = 5;
         elseif currentmcflag == 6
            mcflag = 1;
        elseif currentmcflag == 7
            mcflag = 0;
        end
        H(8) = mcflag;
       
        
        if ~isempty(find(currentmcflag == [1 2 3 4]))
            H(14) =NaN; H(15) = NaN; H(16) = NaN; H(21) = NaN; 
            set([enterxfocus enteryfocus enterzfocus enterwaist], 'String',num2str([]),'Enable','off')
            set(infoLabel,'String',sprintf('light beam is changed to %s.',mcflagtxt{currentmcflag}));

        else
            set([enterxfocus enteryfocus enterzfocus enterwaist],'Enable','on')
            H(14) = H(11);H(15) = H(12); H(16) = H(13)+H(20)*tan(asin(NA/1.36)); H(21)= H(20)/2;
            set(enterxfocus,'String', num2str(H(14)));
            set(enteryfocus,'String', num2str(H(15)));
            set(enterzfocus,'String', num2str(H(16)));
            set(enterwaist,'String', num2str(H(21)));  
            set(infoLabel,'String',sprintf('Light beam is changed to %s. Default values were asigned to xfocus, yfocus, zfocus, and waist',mcflagtxt{currentmcflag}));
            
        end
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
            
    end

%--------------------------------------------------------------------------

%%YL: setcontour_callback

    function setcontour_callback(lightdenmenu,eventdata)
        
        if get(lightdenmenu,'Value') == 1
            defaultcontours = sprintf('%4s %4s %4s %4s','0.1', '1','10','100');
            set(enterlightden,'string',defaultcontours,'enable','off');
            conFLAG = 0;   % use default contours
        elseif get(lightdenmenu,'Value') == 2
            defaultcontours = sprintf('%4s',num2str(conVALUE));
            set(enterlightden,'string',defaultcontours,'enable','on');
            conVALUE = str2num(get(enterlightden,'string'));
            conFLAG = 1;    % use the specified contour
        end
        
    end
%%YL: entercontour_callback

%--------------------------------------------------------------------------
    function entercontour_callback(enterlightden,eventdata)
  
        if get(lightdenmenu,'Value') == 2
            defaultcontours = sprintf('%4s',get(enterlightden,'string'));
            conVALUE = str2num(defaultcontours);
            conFLAG = 1;    % use the specified contour
        end
        
    end


%--------------------------------------------------------------------------
% callback function for simatlas
    function loadATLAS(load3datlas,eventdata)
        VOL2 = permute(VOL,[2 1 3]);
        xe = round(H(11)/dx+1/2+Nx/2);
        ye= round(H(12)/dy+1/2 + Ny/2);
        ze = round(H(13)/dz+1/2);
        %         if ishghandle(sbfigh)
        %             close(sbfigh)
        %         end
        [sbxyz sbfig] = SliceBrowser(VOL2,xe,ye,ze); % Nx, Ny, Nz
        sbfigh = sbfig; %figure(sbfig);
        disp(sprintf('%s',sbxyz))
        %YL:add  a push putton on the atals to check/update the launching position
        enterxyz = uicontrol('Parent',sbfigh,'Style','pushbutton','String','Use current XYZ as launching position','Visible','on',...
            'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized',...
            'Position',[0.468 .066 .398 .054],'callback','ClickedCallback','Callback', {@enterXYZ_Callback});
        
        set(infoLabel,'String','The 3D atlas currently used here is introduced in the paper authored by  Y. Ma, et al in Neuroscience 135, 1203–1215 (2005)');
        
        
    end



%--------------------------------------------------------------------------
% callback function for SIMname
    function textbox_SIM(object,eventdata)
        usr_input = get(enterSIM,'String');
        set(enterSIM,'UserData',usr_input);
        SIMname = usr_input;
        disp(['Simname = ' SIMname])
        [T VOL H SIMname] = loadmc(sprintf('%s',SIMname),librdir);
        set(infoLabel,'String',sprintf('Simulation named %s.mci is loaded.',SIMname));
        set(callsavemci,'backgroundColor','g','String','Click to Save')
        set(callGO,'backgroundColor','g','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

%--------------------------------------------------------------------------
% callback function for POWER
    function textbox_pwr(object,eventdata)
        usr_input = get(enterpwr,'String');
        usr_input = str2double(usr_input);
        set(enterpwr,'UserData',usr_input);
        pwr = usr_input;
        H(23) = pwr;
        set(infoLabel,'String',sprintf('Power from the fiber tip is set to %3.1f mw.',pwr));
%         set(callsavemci,'backgroundColor','r','String','Click to Save')
%         set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
       
    end

%--------------------------------------------------------------------------
% callback function for WAVELENGTH
    function textbox_nm(object,eventdata)
        usr_input = get(enterNM,'String');
        usr_input = str2double(usr_input);
        set(enterNM,'UserData',usr_input);
        nm = usr_input;
        H(22) = nm;
        set(infoLabel,'String',sprintf('Light wavelength is set to %d nm.', nm));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
    end

%--------------------------------------------------------------------------
% callback function for TIME
    function timeorphoton_callback(object, eventdata)
        
        if get(min_numChoice,'Value') == 1   % set time
           timeFLAG = 1 ;
           time_photons = 1;                % set default simulation time to 1 minute
           H(1) = time_photons;  H(27) = timeFLAG;
           set(entertime_photon,'String',num2str(time_photons))
           set(infoLabel,'String',sprintf('Simulation time is set to default %3.2f minute(s).',time_photons));
           set(callsavemci,'backgroundColor','r','String','Click to Save')
        elseif get(min_numChoice,'Value') == 2 % set photons
            timeFLAG = 2 ;
            time_photons = 100000;           % set default photon number to 100,000
            H(1) = time_photons;  H(27) = timeFLAG;
            set(entertime_photon,'String',num2str(time_photons))
            set(infoLabel,'String',sprintf('The number of photons is set to default %d',time_photons));
            set(callsavemci,'backgroundColor','r','String','Click to Save')
        end
    
    end

%--------------------------------------------------------------------------

% callback function for TIME_PHOTON
    function textbox_min(object,eventdata)
        usr_input = get(entertime_photon,'String');
        usr_input = str2double(usr_input);
        set(entertime_photon,'UserData',usr_input);
        time_photons = usr_input;
        H(1) = time_photons;
        if timeFLAG == 1
            set(infoLabel,'String',sprintf('Simulation time is set to %3.2f minute(s).',time_photons));
        elseif timeFLAG == 2
            set(infoLabel,'String',sprintf('The number of photons is set to %d.',time_photons));
        end
                 
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

%--------------------------------------------------------------------------
% callback function for fibe radius
    function textbox_radius(object,eventdata)
        usr_input = get(enterradius,'String');
        usr_input = str2double(usr_input);
        set(enterradius,'UserData',usr_input);
        radius = usr_input;
        H(20) = radius;
        set(infoLabel,'String',sprintf('Fiber core radius is set to %3.1f um.',radius*1e4));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

%--------------------------------------------------------------------------
% callback function for fibe NA
    function textbox_NA(object,eventdata)
        usr_input = get(enterNA,'String');
        usr_input = str2double(usr_input);
        set(enterNA,'UserData',usr_input);
        NA = usr_input;
        H(24) = NA;
        H(25) = asin(NA/1.36);
        set(infoLabel,'String',sprintf('Fiber NA is set to %3.2f.',NA));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

%--------------------------------------------------------------------------
% callback function for xs
    function textbox_xs(object,eventdata)
        usr_input = get(enterxs,'String');
        usr_input = str2double(usr_input);
        set(enterxs,'UserData',usr_input);
        xs = usr_input;
        H(11) = xs;
        set(infoLabel,'String',sprintf('X of the fiber center is set to %6.4f cm.',xs));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end
%--------------------------------------------------------------------------
% callback function for ys
    function textbox_ys(object,eventdata)
        usr_input = get(enterys,'String');
        usr_input = str2double(usr_input);
        set(enterys,'UserData',usr_input);
        ys = usr_input;
        H(12) = ys;
        set(infoLabel,'String',sprintf('Y of the fiber center is set to %6.4f cm.',ys));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end
%--------------------------------------------------------------------------
% callback function for zs
    function textbox_zs(object,eventdata)
        usr_input = get(enterzs,'String');
        usr_input = str2double(usr_input);
        set(enterzs,'UserData',usr_input);
        zs = usr_input;
        H(13) = zs;
        set(infoLabel,'String',sprintf('Z of the fiber center is set to %6.4f cm.',zs));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

%--------------------------------------------------------------------------
  function setmcflag(mcflag)
      
        if     mcflag == 0, currentmcflag = 7;
        elseif mcflag == 1, currentmcflag = 6;
        elseif mcflag == 2, currentmcflag = 4;;
        elseif mcflag == 3, currentmcflag = 3;
        elseif mcflag == 4, currentmcflag = 2;
        elseif mcflag == 5, currentmcflag = 5;
        elseif mcflag == 6, currentmcflag = 1;
        else
            error('Please make sure mcflag is an integer within [0 6]')
        end
        
  end
%--------------------------------------------------------------------------

% callback function for xfocus
    function textbox_xfocus(object,eventdata)
        usr_input = get(enterxfocus,'String');
        usr_input = str2double(usr_input);
        set(enterxfocus,'UserData',usr_input);
        xfocus = usr_input;
        H(14) = xfocus;
        set(infoLabel,'String',sprintf('X of the focus is set to %6.4f cm.',xfocus));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end
%--------------------------------------------------------------------------
% callback function for yfocus
    function textbox_yfocus(object,eventdata)
        usr_input = get(enteryfocus,'String');
        usr_input = str2double(usr_input);
        set(enteryfocus,'UserData',usr_input);
        yfocus = usr_input;
        H(15) = yfocus;
        set(infoLabel,'String',sprintf('Y of the focus is set to %6.4f cm.',yfocus));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end
%--------------------------------------------------------------------------
% callback function for zfocus
    function textbox_zfocus(object,eventdata)
        usr_input = get(enterzfocus,'String');
        usr_input = str2double(usr_input);
        set(enterzfocus,'UserData',usr_input);
        zfocus = usr_input;
        H(16) = zfocus;
        set(infoLabel,'String',sprintf('Z of the focus is set to %6.4f cm.',zfocus));
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end
%--------------------------------------------------------------------------
% callback function for waist at focus
    function textbox_waist(object,eventdata)
        usr_input = get(object,'String');
        usr_input = str2double(usr_input);
        set(object,'UserData',usr_input);
        waist = usr_input;
        H(21) = waist;
        set(infoLabel,'String','Beam waist at focus is set(not valid for de-focused beam).');
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end


%--------------------------------------------------------------------------
%  callback function for enterxyz
    function enterXYZ_Callback(object,eventdata)
        
        pointer3dt = getfield(guidata(sbfig),'pointer3dt');
        disp(sprintf('Previous xs ys zs %s',num2str(H(11:13)')));
        
        usr_input(1) = str2num(sprintf('%5.4f',(pointer3dt(1)-128-1/2)*dx));
        usr_input(2) = str2num(sprintf('%5.4f',(pointer3dt(2)-256-1/2)*dx));
        usr_input(3) = str2num(sprintf('%5.4f',(pointer3dt(3)-1/2)*dx));
        
        
        H(11) = usr_input(1);
        H(12) = usr_input(2);
        H(13) = usr_input(3);
%         H(14) = H(11); % YL: set xfocus = xs
%         H(15) = H(12); % YL: set yfocus = ys
        
        set(enterxs,'String',H(11));
        set(enterys,'String',H(12));
        set(enterzs,'String',H(13));
        
         if ~isempty(find(currentmcflag == [1 2 3 4]))
            set([enterxfocus enteryfocus enterzfocus enterwaist], 'String',num2str([]),'Enable','off')
            set(infoLabel,'String',sprintf('Launching position was updated to [ %s ].',num2str(H(11:13)')));

        else
            set([enterxfocus enteryfocus enterzfocus enterwaist],'Enable','on')
            H(14) = H(11);H(15) = H(12); H(16) = H(13)+H(20)*tan(asin(NA/1.36)); H(21)= H(20)/2;
            set(enterxfocus,'String', num2str(H(14)));
            set(enteryfocus,'String', num2str(H(15)));
            set(enterzfocus,'String', num2str(H(16)));
            set(enterwaist,'String', num2str(H(21)));  
            set(infoLabel,'String',sprintf('Launching position was updated to [ %s ]. \n . Default values were asigned to xfocus, yfocus, zfocus, and waist',num2str(H(11:13)')));
            
        end
      
        set(callsavemci,'backgroundColor','r','String','Click to Save')
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        simFLAG = 0;
        set(callOUTPUT,'Enable','off')
        
    end


%%%%%%
% functions
%%%
%--------------------------------------------------------------------------
% callback function for loadmci
    function loadSIM(callloadsim,eventdata)
        
        %%YL: load different parameters,but save as mc.mci
        [simpName simpPath] = uigetfile({'*.mci';'*.*'},...
            'Load parameters via .mci file','./','MultiSelect','off');
         
        flagF = (exist(strrep(simpName,'.mci','_F.bin'))>0);
        if flagF
            copyfile(fullfile(simpPath,simpName),fullfile(librdir,'mc.mci'));
            copyfile(fullfile(simpPath,strrep(simpName,'.mci','_F.bin')),fullfile(librdir,'mc_F.bin'));
            disp(['Copying ' strrep(simpName,'.mci','_F.bin')])

        else
            copyfile(fullfile(simpPath,simpName),fullfile(librdir,'mc.mci'));
            set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
            set(callOUTPUT,'Enable','off')
            set(infoLabel,'String','No F.bin file exists yet.');
            disp('No F.bin file exists yet.')
            
        end
        [~,~,H,~] = loadmc('mc',librdir); % <-------- create .mci file -------
        if H(16) ==  1.0000e+12
            H(16) = inf;
        end
        
        set(infoLabel,'String',sprintf('%s was loaded and mc.mci was updated',fullfile(simpPath,simpName)));
      % <-------- create .mci file -------
        %         set(infoLabel,'String','mc.mci file saved.');
        time_photons = H(1); % time in simulation
        Nx      = H(2);
        Ny      = H(3);
        Nz      = H(4);
        dx      = H(5);
        dy      = H(6);
        dz      = H(7);
        mcflag  = H(8);
        launchflag = H(9);
        boundaryflag = H(10);
        xs      = H(11);
        ys      = H(12);
        zs      = H(13);
        xfocus  = H(14);
        yfocus  = H(15);
        zfocus  = H(16);
        ux0     = H(17);
        uy0     = H(18);
        uz0     = H(19);
        radius  = H(20);
        waist   = H(21);
        nm      = H(22);
        pwr     = H(23);
     %%YL: add NA and thmax  
        NA      = H(24);
        thmax   = H(25);
        
        Nt      = H(26);
        timeFLAG = H(27); %YL
        set(enterSIM,'string',simpName);
        set(enterpwr,'string',num2str(pwr));
        set(enterNM,'string',num2str(nm));
        set(entertime_photon,'string',num2str(time_photons));
        set(enterxs,'string',num2str(xs));
        set(enterys,'string',num2str(ys));
        set(enterzs,'string',num2str(zs));
        setmcflag(mcflag); set(beamchoice,'Value',currentmcflag);
        set(enterxfocus,'string',num2str(xfocus));
        set(enteryfocus,'string',num2str(yfocus));
        set(enterzfocus,'string',num2str(zfocus));
        set(enterradius,'string',num2str(radius));
        set(enterwaist,'string',num2str(waist));
        
        if ~isempty(find(currentmcflag == [1 2 3 4]))
            set([enterxfocus enteryfocus enterzfocus enterwaist], 'String',num2str([]),'Enable','off')
        else
            set([enterxfocus enteryfocus enterzfocus enterwaist],'Enable','on')
        end
        
       if flagF == 0
           set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
           set(callOUTPUT,'Enable','off')
       else
           set(callGO,'backgroundColor','g','Enable','on')
           set(callOUTPUT,'backgroundColor','g','Enable','on')
           
       end
        
    end

%--------------------------------------------------------------------------
% callback function for savemci
    function savemci(callsavemci,eventdata)
        if H(16) ==  inf
            H(16) = 1.0000e+12;
        end
        saveHmci('mc', librdir, H,nm, tissue,originalname); % <-------- create .mci file -------
        set(callsavemci,'backgroundColor','g','String','Saved')
        c = round(clock);
        switch COMP
            case 1 % 'Mac'
                bkupname = sprintf('mc_%d-%d_%d_%d_%d',c(4),c(5),c(2),c(3),c(1))% SLJ
            case 0 % 'Win'
                bkupname = sprintf('mc_%d-%d_%d_%d_%d',c(4),c(5),c(2),c(3),c(1)) %
        end
        %
        
        saveHmci(bkupname, librdir, H,nm, tissue,originalname); % <-------- create .mci file -------
        set(infoLabel,'String','mc.mci file saved.');
        set(enterSIM,'string',bkupname)
        if H(16) == 1.0000e+12
            H(16) = inf;
        end
    end

%--------------------------------------------------------------------------
% callback function for GO button
    function GO(callGO,eventdata)
        if flagOUTGO
            
            if exist(fullfile(librdir,'time_min.csv'),'file');
               delete(fullfile(librdir,'time_min.csv'),'file');  % delete the previous time file if exists
            end
        
            if H(16) ==  inf
                H(16) = 1.0000e+12;
            end
            saveHmci('mc', librdir, H,nm, tissue,originalname); % <-------- create .mci file -------
            set(callGO,'backgroundColor','g')
            set(callsavemci,'backgroundColor','g','String','Saved')
            set(infoLabel,'String',sprintf('mc.mci file saved.\nSimulation is ongoing. Check the progress in the commond window.'));
            pause(.1)
            disp('gomcxyzOGS mc')
            set([enterpwr enterSIM enterNM min_numChoice entertime_photon enterradius enterNA enterxs enterys enterzs],'Enable','off')
            set([beamchoice enterxfocus enteryfocus enterzfocus enterwaist enterxyz],'Enable','off')
            set([callsavemci loadsimulation lightdenmenu enterlightden],'Enable','off')

            switch COMP
                case 1 % 'Mac'
%                     status = system('./gomcxyz OGSout/mc') % SLJ
                   
                    pdtemp = pwd; disp(sprintf('Before simulation, folder is at %s',pdtemp))
                    cd(librdir)
                    
                    system('./gomcxyzOGS mc') % SLJ
                        
                    cd(pdtemp);disp(sprintf('After simulation, folder is at %s',pdtemp))
                    
                    set(infoLabel,'String',' mcxyz is done');

                case 0 % 'Win'
%                   status = system('gomcxyzMay8 OGSout\mc &') % SLJ
                    pdtemp = pwd; disp(sprintf('Before simulation, folder is at %s',pdtemp))
                    cd(librdir)
                    system('gomcxyzOGS mc &'); % SLJ
                    set(infoLabel,'String',' Simultion  is ongoing.');
                   
                    if timeFLAG == 1
                        time_min = time_photons;
                        pause(round(60*time_min)+2)
                        set(infoLabel,'String',' mcxyz is done');
                        cd(pdtemp); disp(sprintf('After simulation, folder is at %s',pdtemp))
                    elseif timeFLAG == 2
                         Nphotons = time_photons;
                         pause(10)   % wait for 10 seconds to make sure the csv file is created 
                       if exist(fullfile(librdir,'time_min.csv'),'file');
                           time_min = importdata(fullfile(librdir,'time_min.csv'));
                           set(infoLabel,'String',sprintf(' Simultion  is in progress. Estimated time = %4.3f min for %d photons',time_min,Nphotons));
                           pause(round(60*time_min))
                           cd(pdtemp); disp(sprintf('After simulation, folder is at %s',pdtemp))
                          set(infoLabel,'String',' mcxyz is done');
                        
                       else
                        set(infoLabel,'String',' Simultion  is ongoing,check the commond window to make sure the simulation is done before proceed.');
                        cd(pdtemp); disp(sprintf('After simulation, folder is at %s',pdtemp))
                       end
                    end
                   

            end
     
            
            set(callGO,'Enable','on')
            set([enterpwr lightdenmenu enterlightden],'Enable','on')
           
            if H(16) == 1.0000e+12
                H(16) = inf;
            end
            
            
            set(callOUTPUT,'Enable','on')
                           
                       
        else
            disp('GO button disabled while editing figure.')
        end
        
    end

%--------------------------------------------------------------------------
% callback function for output button
    function OUTPUT(callOUTPUT,eventdata)
        if flagOUTGO
            flagOUTGO = 0;
            grayclr = [1 1 1]*.4;
            set(callGO,'backgroundColor',grayclr)
            set(callOUTPUT,'backgroundColor',grayclr)
            
            if H(16) ==  inf
                H(16) = 1.0000e+12;
            end
            %             saveHmci('mc', './', H,nm, tissue,originalname); % <-------- create .mci file -------
            saveHmci('mc', librdir, H,nm, tissue,originalname); % YL: add sub-folder name 'OGSout' <-------- create .mci file -------
            
            set(callsavemci,'backgroundColor','g','String','Saved')
            set(infoLabel,'String','mc.mci saved');
            pause(0.3)
            %             flagF = exist('mc_F.bin');
            flagF = exist(fullfile(librdir,'mc_F.bin'));
            if flagF==0
                set(infoLabel,'String','No mc_F.bin file exists yet.');
                disp('No mc_F.bin file exists yet.')
            else
                %% % Load Fluence rate F(y,x,z)
                set([enterpwr lightdenmenu enterlightden],'Enable','off')

                set(infoLabel,'String','loading mc_F.bin');
                mcFname = fullfile(librdir,'mc_F.bin');
                disp(['loading ' mcFname])
                tic
                fid = fopen(mcFname, 'rb');
                [Data count] = fread(fid, Ny*Nx*Nz, 'float');
                fclose(fid);
                toc
                F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
                ss =            'to select x,y,z   --> click on figure then hit RETURN.';
                ss = strvcat(ss,'to select view    --> click on box,then hit RETURN.');
                ss = strvcat(ss,'to RETURN to here --> click outside image then hit RETURN.');
                set(infoLabel,'String',ss);
                bj = 3;     % choose to view fluence contours
                
                %% YL: add the absolute threshold for the desired light power density
                
                
                [ix,iy,iz] = lookOGSoutput(F,VOL,T,bj,conFLAG,conVALUE,librdir);
                
                set(infoLabel,'String',note1)
                set([enterpwr lightdenmenu enterlightden],'Enable','on')
            end
            flagOUTGO = 1;
            set(callGO,'backgroundColor','g')
            set(callOUTPUT,'backgroundColor','g')
        else
            disp('OUTPUT button disabled while editing figure.')
        end
    end

%--------------------------------------------------------------------------
% callback function for CHECKMCI
    function checkmci(callcheckmci,eventdata)
        home
%         reportHmci('mc',1,librdir);
        reportHmci(get(enterSIM,'String'),1,librdir);  % YL1024
        set(infoLabel,'String','Parameters for this simulation are listed in the command window.');
    end

%--------------------------------------------------------------------------
% callback function for DEBUG
    function debugON(calldebugON,eventdata)
        fprintf('debug on\n')
        set(infoLabel,'String',sprintf('debug on.\nType "dbquit" in command window to resume.\nHit <Reset> to clear this message.'))
        keyboard
    end

%--------------------------------------------------------------------------
% callback function for RESET
   function resetOGS(~,~)
       try 
           close(figure(501))
       catch
           disp('Output window does not exist while resetting OptogenSIM')
       end
       
       try
       
            OptogenSIM
       catch 
            cd ..
            disp('Simulation may be interrupted, return to the original directory to start')
            OptogenSIM
       end
%         eval([optogensim])
    end

% disp('ready')
end % ends optogensim function
