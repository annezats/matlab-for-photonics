function varargout = light_polarization(varargin)
% LIGHT_POLARIZATION MATLAB code for light_polarization.fig
%      LIGHT_POLARIZATION, by itself, creates a new LIGHT_POLARIZATION or raises the existing
%      singleton*.
%
%      H = LIGHT_POLARIZATION returns the handle to a new LIGHT_POLARIZATION or the handle to
%      the existing singleton*.
%
%      LIGHT_POLARIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIGHT_POLARIZATION.M with the given input arguments.
%
%      LIGHT_POLARIZATION('Property','Value',...) creates a new LIGHT_POLARIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before light_polarization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to light_polarization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help light_polarization

% Last Modified by GUIDE v2.5 06-Jul-2020 16:18:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @light_polarization_OpeningFcn, ...
                   'gui_OutputFcn',  @light_polarization_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before light_polarization is made visible.
function light_polarization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to light_polarization (see VARARGIN)

% Choose default command line output for light_polarization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes light_polarization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = light_polarization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function amplitude1_Callback(hObject, eventdata, handles)
% hObject    handle to amplitude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amplitude1 as text
%        str2double(get(hObject,'String')) returns contents of amplitude1 as a double
amp1 = str2num(get(handles.amplitude1,'String'));
setappdata(0,'amplitude1',amp1);


% --- Executes during object creation, after setting all properties.
function amplitude1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amplitude2_Callback(hObject, eventdata, handles)
% hObject    handle to amplitude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amplitude2 as text
%        str2double(get(hObject,'String')) returns contents of amplitude2 as a double
amp2 = str2num(get(handles.amplitude2,'String'));
setappdata(0,'amplitude2',amp2);




% --- Executes during object creation, after setting all properties.
function amplitude2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function omega1_Callback(hObject, eventdata, handles)
% hObject    handle to omega1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of omega1 as text
%        str2double(get(hObject,'String')) returns contents of omega1 as a double
omega_1 = str2num(get(handles.omega1,'String'));
setappdata(0,'omega1',omega_1);


% --- Executes during object creation, after setting all properties.
function omega1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function omega2_Callback(hObject, eventdata, handles)
% hObject    handle to omega2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of omega2 as text
%        str2double(get(hObject,'String')) returns contents of omega2 as a double
omega_2 = str2num(get(handles.omega2,'String'));
setappdata(0,'omega2',omega_2);




% --- Executes during object creation, after setting all properties.
function omega2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to omega2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


cla
% constants
k= 2*pi; % frequency
z=0:0.01:3; % direction of propogation

A1 = getappdata(0,'amplitude1'); %Amplitudes
A2 = getappdata(0,'amplitude2');

O1 = getappdata(0,'omega1')*2* pi; %Omegas
O2 = getappdata(0,'omega2')* 2*pi;

Ex= A1 * sin(k*z-O1); % Electric Fields
Ey= A2 * sin(k*z-O2);

omega = mod(abs(O1 - O2),pi);

amplitude = A1 - A2;

circle = pi/2;


    
hold on
%----total polarization
plot3(z,Ex,Ey,'LineWidth',2)
%----x-y projection
plot3(z*0 + 3,Ex,Ey,'LineWidth',2)
%-----x and y components
plot3(z,Ex,Ey*0,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,Ex,Ey*0,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
plot3(z,Ex*0,Ey,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,Ex*0,Ey,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
%----making a line through the origin
plot(z,z*0,'LineWidth',2)
%---Display
grid on
xlim([0 3]);

if A1 > A2
    ylim([-1*A1 A1]);
    zlim([-1*A1 A1]);
else
    ylim([-1*A2 A2]);
    zlim([-1*A2 A2]);
end
    
xlabel('Z/wavelength'); 
ylabel('Ex');
zlabel('Ey');
view(-40,40);
hold off

if omega == 0  %Linear Polarization
    title('\fontsize{18} Linear Polarization')
   
else
    if omega == circle && amplitude == 0 %Circular Polarization
        title('\fontsize{18} Circular Polarization')

    else %Elliptical Polarization
        title('\fontsize{18} Elliptical Polarization')
        
    end
end

       


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
% Getting Values


% Funded by the CUNY Summer Undergrad Research Program
% https://www.cuny.edu/research/student-resources/for-students/cuny-summer-undergraduate-research-program/
% Created by : Natascha Krishnanand, Anne Zats and Anton Kyrylenko