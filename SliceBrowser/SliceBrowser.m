% ======================================================================
%> SLICEBROWSER M-file for SliceBrowser.fig
%>       SliceBrowser is an interactive viewer of 3D volumes, 
%>       it shows 3 perpendicular slices (XY, YZ, ZX) with 3D pointer.
%>   Input:  a) VOLUME - a 3D matrix with volume data
%>           b) VOLUME - a 4D matrix with volume data over time
%>   Control:
%>       - Clicking into the window changes the location of 3D pointer.
%>       - 3D pointer can be moved also by keyboard arrows.
%>       - Pressing +/- will switch to next/previous volume.
%>       - Pressing 1,2,3 will change the focus of current axis.
%>       - Pressing 'e' will print the location of 3D pointer.
%>       - Pressing 'c' switches between color-mode and grayscale.
%>       - Pressing 'q' switches scaling of axis (equal/normal).
%>   Example of usage:
%>       load mri.dat
%>       volume = squeeze(D);
%>       SliceBrowser(volume);
%>
%> Author: Marian Uhercik, CMP, CTU in Prague
%> Web: http://cmp.felk.cvut.cz/~uhercik/3DSliceViewer/3DSliceViewer.htm
%> Last Modified by 21-Jul-2011
% ======================================================================

function varargout = SliceBrowser(varargin)

% Documentation generated GUIDE:
%
%SLICEBROWSER M-file for SliceBrowser.fig
%      SLICEBROWSER, by itself, creates a new SLICEBROWSER or raises the existing
%      singleton*.
%
%      H = SLICEBROWSER returns the handle to a new SLICEBROWSER or the handle to
%      the existing singleton*.
%
%      SLICEBROWSER('Property','Value',...) creates a new SLICEBROWSER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SliceBrowser_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SLICEBROWSER('CALLBACK') and SLICEBROWSER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SLICEBROWSER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SliceBrowser

% Last Modified by GUIDE v2.5 12-May-2015 19:36:02


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SliceBrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @SliceBrowser_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% % YL: handles structure
% SliceBrowserFigure: 1
% pointer3d_info: 716.0291
% Subplot3: 711.0291
% Subplot2: 706.0291
% Subplot1: 701.0291
% output: 1
% volume: [256x512x256 double]
% axis_equal: 1
% color_mode: 0
% pointer3dt: [150 257 64 1]
% vol_sz: [256 512 256 1]
% last_axis_id: 0
% 

% --- Executes just before SliceBrowser is made visible.
function SliceBrowser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SliceBrowser
handles.output = hObject;

% handles.SliceBrowserFigure = varargin{5};  % YL
% gcf110 = varargin{5}; % YL

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SliceBrowser wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if (length(varargin) <=0)
    error('Input volume has not been specified.');
end;
volume = varargin{1};
if (ndims(volume) ~= 3 && ndims(volume) ~= 4)
    error('Input volume must have 3 or 4 dimensions.');
end;
handles.volume = volume;

 %handles.axis_equal = 0;
handles.axis_equal = 1; %YL

handles.color_mode = 1;
if (size(volume,4) ~= 3)
    handles.color_mode = 0;
end;

% set main wnd title
% set(gcf, 'Name', 'Slice Viewer: 3D atlas')
%YL
sz0 = get(0,'screensize'); 
sw0 = sz0(3);  
sh0 = sz0(4);
whr = sw0/sh0; % screen width / height ratio
pg = [0.01, 0.25 0.25 0.68];   % position of control gui
% pg = [0.01, 0.25 0.25 0.68];   % position of control gui


% gcf110 = gcf;%figure(110);
set(hObject,'Resize','on','Units','pixels','Position',[sw0*(pg(1)+pg(3)+0.012) sh0*(pg(2)+pg(4)/2*0) sh0*pg(4)*2.2*0.50 sh0*pg(4)*2.2*0.50],'Visible','on',...
    'MenuBar','none','name','3D mouse brain atlas','NumberTitle','off','UserData',0);      % enable the Menu bar so that to explore the intensity value
   
% init 3D pointer
vol_sz = size(volume); 
if (ndims(volume) == 3)
    vol_sz(4) = 1;
end;
pointer3dt = floor(vol_sz/2)+1;

%YL: use the input from the main GUI
pointer3dt(1)= varargin{2};
pointer3dt(2)= varargin{3};
pointer3dt(3)= varargin{4};

handles.pointer3dt = pointer3dt;
% handles.pointer3dt

handles.vol_sz = vol_sz;

% YL: add three text boxes to specify the X, Y, Z position of  the
% launching poisiton
% xyzlabel = uicontrol('Parent',gcf,'Style','text','String',['Current:   X   Y   Z'],...
%     'FontUnits','normalized','FontSize',.3,'Units','normalized','Position',[.595 .15 .3 .1]);
% YL: use the input from the main GUI
dx = 0.0047;
xe = varargin{2};
ye = varargin{3};
ze = varargin{4}; 
xsc = sprintf('%5.4f',(handles.pointer3dt(1)-128-1/2)*dx);
ysc = sprintf('%5.4f',(handles.pointer3dt(2)-256-1/2)*dx);
zsc = sprintf('%5.4f',(handles.pointer3dt(3)-1/2)*dx);

set(handles.getXYZinvox,'String',sprintf('%d     %d     %d',xe,ye,ze));
set(handles.getXYZincm ,'String', sprintf('%s     %s     %s', xsc, ysc, zsc));
plot3slices(hObject, handles);
% stores ID of last axis window
% (0 means that no axis was clicked yet)
handles.last_axis_id = 0;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = SliceBrowser_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output; %YL
varargout{1} = get(handles.getXYZincm,'String');
varargout{2} = handles.SliceBrowserFigure;
%varargout{2} = get(handles.SliceBrowserFigure,'Number');


% varargout{2} = handles.pointer3dt;

% --- Executes on mouse press over axes background.
function Subplot3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XY slice
% This object contains the ZX slice  %YL

% disp('Subplot3:BtnDown');
pt=get(hObject,'currentpoint');
% xpos=round(pt(1,2)); ypos=round(pt(1,1));
xpos=round(pt(1,1)); ypos=round(pt(1,2));  % YL

zpos = handles.pointer3dt(3);
tpos = handles.pointer3dt(4);
handles.pointer3dt = [xpos ypos zpos tpos];
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
plot3slices(hObject, handles);

% % store this axis as last clicked region
 handles.last_axis_id = 3;
% % Update handles structure
 guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function Subplot1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the YZ slice
% This object contains the ZY slice % YL

% disp('Subplot1:BtnDown');
pt=get(hObject,'currentpoint');
% xpos=round(pt(1,2)); zpos=round(pt(1,1));
xpos=round(pt(1,1)); zpos=round(pt(1,2)); % YL

ypos = handles.pointer3dt(2);
tpos = handles.pointer3dt(4);
handles.pointer3dt = [xpos ypos zpos tpos];
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
plot3slices(hObject, handles);

% store this axis as last clicked region
handles.last_axis_id = 1;
% Update handles structure
guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function Subplot2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Subplot3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This object contains the XZ slice
% This object contains the YX slice  % YL

% disp('Subplot2:BtnDown');
pt=get(hObject,'currentpoint');
% zpos=round(pt(1,2)); ypos=round(pt(1,1));
% xpos = handles.pointer3dt(1);
zpos=round(pt(1,2)); ypos=round(pt(1,1)); %YL
xpos = handles.pointer3dt(1);


tpos = handles.pointer3dt(4);
handles.pointer3dt = [xpos ypos zpos tpos];
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
plot3slices(hObject, handles);

% store this axis as last clicked region
handles.last_axis_id = 2;
% Update handles structure
guidata(hObject, handles);

% --- Executes on key press with focus on SliceBrowserFigure and no controls selected.
function SliceBrowserFigure_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to SliceBrowserFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disp('SliceBrowserFigure_KeyPressFcn'); %YL
% curr_char = int8(get(gcf,'CurrentCharacter'));
curr_char = int8(get(hObject,'CurrentCharacter')); %YL

if isempty(curr_char)
    return;
end;


xpos = handles.pointer3dt(1);
ypos = handles.pointer3dt(2);
zpos = handles.pointer3dt(3); 
tpos = handles.pointer3dt(4); 
% Keys:
% - up:   30
% - down:   31
% - left:   28
% - right:   29
% - '1': 49
% - '2': 50
% - '3': 51
% - 'e': 101
% - plus:  43
% - minus:  45
switch curr_char
    case 99 % 'c'
        handles.color_mode = 1 - handles.color_mode;
        if (handles.color_mode ==1 && size(handles.volume,4) ~= 3)
            handles.color_mode = 0;
        end;
        
    case 113 % 'q'
        handles.axis_equal = 1 - handles.axis_equal;
        
    case 30
%         switch handles.last_axis_id
%             case 1
%                 xpos = xpos -1;
%             case 2
%                 xpos = xpos -1;
%             case 3
%                 zpos = zpos -1;
%             case 0
%         end;
          switch handles.last_axis_id
            case 1
                zpos = zpos -1;
            case 2
                zpos = zpos -1;
            case 3
                ypos = ypos -1;
            case 0
        end;
    case 31
%         switch handles.last_axis_id
%             case 1
%                 xpos = xpos +1;
%             case 2
%                 xpos = xpos +1;
%             case 3
%                 zpos = zpos +1;
%             case 0
%         end;
          switch handles.last_axis_id
            case 1
                zpos = zpos +1;
            case 2
                zpos = zpos +1;
            case 3
                ypos = ypos +1;
            case 0
        end;
    case 28
%         switch handles.last_axis_id
%             case 1
%                 ypos = ypos -1;
%             case 2
%                 zpos = zpos -1;
%             case 3
%                 ypos = ypos -1;
%             case 0
%         end;
          switch handles.last_axis_id
            case 1
                xpos = xpos -1;
            case 2
                ypos = ypos -1;
            case 3
                xpos = xpos -1;
            case 0
        end;
    case 29
%         switch handles.last_axis_id
%             case 1
%                 ypos = ypos +1;
%             case 2
%                 zpos = zpos +1;
%             case 3
%                 ypos = ypos +1;
%             case 0
%         end;
       switch handles.last_axis_id
            case 1
                xpos = xpos +1;
            case 2
                ypos = ypos +1;
            case 3
                xpos = xpos +1;
            case 0
        end;
    case 43
        % plus key
        tpos = tpos+1;
    case 45
        % minus key
        tpos = tpos-1;
    case 49
        % key 1
        handles.last_axis_id = 1;
    case 50
        % key 2
        handles.last_axis_id = 2;
    case 51
        % key 3
        handles.last_axis_id = 3;
    case 101
        disp(['[' num2str(xpos) ' ' num2str(ypos) ' ' num2str(zpos) ' ' num2str(tpos) ']']);
    otherwise
        return
end;
handles.pointer3dt = [xpos ypos zpos tpos];
handles.pointer3dt = clipointer3d(handles.pointer3dt,handles.vol_sz);
plot3slices(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Plots all 3 slices XY, YZ, XZ into 3 subplots
function [sp1,sp2,sp3] = plot3slices(hObject, handles)
% %YL 
% pointer3d     3D coordinates in volume matrix (integers)
xe = handles.pointer3dt(1);
ye = handles.pointer3dt(2);
ze = handles.pointer3dt(3); 
dx = 0.0047;
xsc = sprintf('%5.4f',(handles.pointer3dt(1)-128-1/2)*dx);
ysc = sprintf('%5.4f',(handles.pointer3dt(2)-256-1/2)*dx);
zsc = sprintf('%5.4f',(handles.pointer3dt(3)-1/2)*dx);

% disp(sprintf('current position: x = %d, y = %d, z = %d',xe,ye,ze));
% handles.pointer3dt;
% size(handles.volume);
value3dt = handles.volume(handles.pointer3dt(1), handles.pointer3dt(2), handles.pointer3dt(3), handles.pointer3dt(4));

set(handles.getXYZincm ,'String', sprintf('%s %s %s', xsc, ysc, zsc));
set(handles.getXYZinvox,'String',sprintf('%d     %d     %d',xe,ye,ze));
set(handles.points3d_intensity,'String',num2str(value3dt));

guidata(hObject, handles);
% 
if (handles.color_mode ==1)
%     sliceXY = squeeze(handles.volume(:,:,handles.pointer3dt(3),:));
    sliceYZ = squeeze(handles.volume(handles.pointer3dt(1),:,:,:));
%     sliceXZ = squeeze(handles.volume(:,handles.pointer3dt(2),:,:));
    sliceYX = (squeeze(handles.volume(:,:,handles.pointer3dt(3),:)))'; % YX 
    sliceZX = (squeeze(handles.volume(:,handles.pointer3dt(2),:,:)))'; %YL


%     max_xyz = max([ max(sliceXY(:)) max(sliceYZ(:)) max(sliceXZ(:)) ]);
%     min_xyz = min([ min(sliceXY(:)) min(sliceYZ(:)) min(sliceXZ(:)) ]);
      max_xyz = max([ max(sliceYX(:)) max(sliceYZ(:)) max(sliceZX(:)) ]); % YL
      min_xyz = min([ min(sliceYX(:)) min(sliceYZ(:)) min(sliceZX(:)) ]);  % YL
    clims = [ min_xyz max_xyz ];
else
%     sliceXY = squeeze(handles.volume(:,:,handles.pointer3dt(3),handles.pointer3dt(4)));
    sliceYZ = squeeze(handles.volume(handles.pointer3dt(1),:,:,handles.pointer3dt(4)));
%     sliceXZ = squeeze(handles.volume(:,handles.pointer3dt(2),:,handles.pointer3dt(4)));
     sliceYX = (squeeze(handles.volume(:,:,handles.pointer3dt(3),handles.pointer3dt(4))))';
     sliceZX = (squeeze(handles.volume(:,handles.pointer3dt(2),:,handles.pointer3dt(4))))';

%     max_xyz = max([ max(sliceXY(:)) max(sliceYZ(:)) max(sliceXZ(:)) ]);
%     min_xyz = min([ min(sliceXY(:)) min(sliceYZ(:)) min(sliceXZ(:)) ]);
    max_xyz = max([ max(sliceYX(:)) max(sliceYZ(:)) max(sliceZX(:)) ]); %YL
    min_xyz = min([ min(sliceYX(:)) min(sliceYZ(:)) min(sliceZX(:)) ]); %YL

    clims = [ min_xyz max_xyz ];
end;
sliceZY = squeeze(permute(sliceYZ, [2 1 3]));

%YL:subplot doesnot work for custom positions
% sp3 = subplot(2,2,3); % YL: change (2,2,1) to (2,2,3) 
axes(handles.Subplot3);
cla;

%colorbar;
% imagesc(sliceXY, clims);
% title('Slice XY');
% ylabel('X');xlabel('Y');
% line([handles.pointer3dt(2) handles.pointer3dt(2)], [0 size(handles.volume,1)]);
% line([0 size(handles.volume,2)], [handles.pointer3dt(1) handles.pointer3dt(1)]);
imagesc(sliceYX, clims); 
% title('Slice YX');
% ylabel('Y');xlabel('X');
line([handles.pointer3dt(1) handles.pointer3dt(1)], [0 size(handles.volume,2)]);
line([0 size(handles.volume,1)], [handles.pointer3dt(2) handles.pointer3dt(2)]);


%set(allchild(gca),'ButtonDownFcn',@Subplot1_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','SliceBrowser(''Subplot3_ButtonDownFcn'',gca,[],guidata(gcbo))');
if (handles.axis_equal == 1)
    axis image; % YL: change from 'axis image' to 'axis image equal'
else
    axis normal;
end;
text(1,260,'Y','fontsize',12,'color','w')
text(127,501,'X','fontsize',12,'color','w')

%YL:subplot doesnot work for custom positions
% sp1 = subplot(2,2,1);  % YL: change (2,2,2) to (2,2,1)
axes(handles.Subplot1);
cla;

% YL
% imagesc(sliceXZ, clims);
% title('Slice XZ');
% ylabel('X');xlabel('Z');

% line([handles.pointer3dt(3) handles.pointer3dt(3)], [0 size(handles.volume,1)]);
% line([0 size(handles.volume,3)], [handles.pointer3dt(1) handles.pointer3dt(1)]);

imagesc(sliceZX, clims);
% title('Slice ZX');
% ylabel('Z');xlabel('X');
line([handles.pointer3dt(1) handles.pointer3dt(1)], [0 size(handles.volume,3)]);
line([0 size(handles.volume,1)], [handles.pointer3dt(3) handles.pointer3dt(3)]);


%set(allchild(gca),'ButtonDownFcn',@Subplot2_ButtonDownFcn);
set(allchild(gca),'ButtonDownFcn','SliceBrowser(''Subplot1_ButtonDownFcn'',gca,[],guidata(gcbo))');
if (handles.axis_equal == 1)
    axis image; % YL: change from 'axis image' to 'axis image equal'
else
    axis normal;
end;
text(127,245,'X','fontsize',12,'color','w')
text(1,127,'Z','fontsize',12,'color','w')


%YL:subplot doesnot work for custom positions
% sp2 = subplot(2,2,2); axis off;% YL: change (2,2,3) to (2,2,2)
axes(handles.Subplot2);
cla;

imagesc(sliceZY, clims); 
% title('Slice ZY');
% ylabel('Z');xlabel('Y');%axis off
line([0 size(handles.volume,2)], [handles.pointer3dt(3) handles.pointer3dt(3)]);
line([handles.pointer3dt(2) handles.pointer3dt(2)], [0 size(handles.volume,3)]);
%set(allchild(gca),'ButtonDownFcn',@Subplot3_ButtonDownFcn);
% set(allchild(gca),'ButtonDownFcn','SliceBrowser(''Subplot3_ButtonDownFcn'',gca,[],guidata(gcbo))');
set(allchild(gca),'ButtonDownFcn','SliceBrowser(''Subplot2_ButtonDownFcn'',gca,[],guidata(gcbo))');  %YL

if (handles.axis_equal == 1)
    axis image; % YL: change from 'axis image' to 'axis image equal'
else
    axis normal;
 end;
text(255,245,'Y','fontsize',12,'color','w')
text(1,127,'Z','fontsize',12,'color','w')

% set(sp1,'position',[0.08 0.69 0.30 0.30]);  %% YL: add position control for subplot1
% set(sp2,'position',[0.39 0.69 0.60 0.30]);  %%YL: add position control for subplot2
% set(sp3,'position',[0.08 0.08 0.30 0.60]);  %% YL: add position control for subplot3, 
% axis([sp1 sp2 sp3],'off');

%  guidata(hObject, handles);



function pointer3d_out = clipointer3d(pointer3d_in,vol_size)
pointer3d_out = pointer3d_in;
for p_id=1:4
    if (pointer3d_in(p_id) > vol_size(p_id))
        pointer3d_out(p_id) = vol_size(p_id);
    end;
    if (pointer3d_in(p_id) < 1)
        pointer3d_out(p_id) = 1;
    end;
end;


function getXYZincm_Callback(hObject,eventdata, handles)
% hObject    handle to getXYZincm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of getXYZincm as text
%        str2double(get(hObject,'String')) returns contents of getXYZincm as a double
%     disp(handles.pointer3dt);
    disp(sprintf('xyz, entered is :%s',get(hObject,'string')));
    usr_input = get(hObject,'String');
    usr_input = str2num(usr_input);
    set(handles.getXYZincm,'UserData',usr_input);
    xsc = usr_input(1);
    ysc = usr_input(2);
    zsc = usr_input(3);
    dx = 0.0047; dy = dx; dz = dx;
    Nx = 256; Ny = 512; Nz = 256;
    xe = round(xsc/dx+1/2+Nx/2);
    ye= round(ysc/dy+1/2 + Ny/2);
    ze = round(zsc/dz+1/2);
    handles.pointer3dt(1:3) = [xe ye ze] ;
    plot3slices(hObject, handles);
