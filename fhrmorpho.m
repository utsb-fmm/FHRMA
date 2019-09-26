% An interface to select a method and display results of various
% Morphological analysis method.
%
%
% USAGE
%    fhrmorpho launch the interface

%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     FHR Morphological Analysis Toolbox  Copyright (C) 2018 Samuel Boudet, Faculté de Médecine et Maïeutique,
%     samuel.boudet@gmail.com
%
%     This file is part of FHR Morphological Analysis Toolbox 
%
%     FHR Morphological Analysis Toolbox  is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     FHR Morphological Analysis Toolbox  is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function varargout = fhrmorpho(varargin)
% FHRMORPHO MATLAB code for fhrmorpho.fig
%      FHRMORPHO, by itself, creates a new FHRMORPHO or raises the existing
%      singleton*.
%
%      H = FHRMORPHO returns the handle to a new FHRMORPHO or the handle to
%      the existing singleton*.
%
%      FHRMORPHO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FHRMORPHO.M with the given input arguments.
%
%      FHRMORPHO('Property','Value',...) creates a new FHRMORPHO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fhrmorpho_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fhrmorpho_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fhrmorpho

% Last Modified by GUIDE v2.5 01-Oct-2018 16:23:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fhrmorpho_OpeningFcn, ...
                   'gui_OutputFcn',  @fhrmorpho_OutputFcn, ...
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


% --- Executes just before fhrmorpho is made visible.
function fhrmorpho_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fhrmorpho (see VARARGIN)

% Choose default command line output for fhrmorpho

handles.listmethod={
    'Ayres*','aamayres(FHR)';
    'Cazares','aamcazares(FHR)';
    'Cazares*','aamcazares(FHR,1)';
    'Houze*','aamhouze(FHR,TOCO,1)';
    'Jimenez','aamjimenez(FHR)';
    'Jimenez*','aamjimenez(FHR,1)';
    'Lu*','aamlu(FHR)';
    'Maeda*','aammaeda(FHR)';
    'Mantel','aammantel(FHR,FHR0)';
    'Mantel*','aammantel(FHR,FHR0,1)';
    'Mongelli*','aammongelli(FHR)';
    'Pardey*','aampardey(FHR)';
    'Taylor','aamtaylor(FHR)';
    'Taylor*','aamtaylor(FHR,1)';
    'Wrobel*','aamwrobel(FHR)';
    'WMFB','aamwmfb(FHR)';
    'Expert consensus','expert';
    };

for i=1:length(handles.listmethod)
    li=floor((i-1)/4);
    col=rem((i-1),4);
    handles.uiMethodToggle(i)=uicontrol(handles.uipanel1,'style','togglebutton',...
        'String',handles.listmethod{i,1},'Position',[10+col*110 140-li*30 100 30]);
end
    
handles.output = hObject;

handles.listPlotter=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fhrmorpho wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fhrmorpho_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,d]=uigetfile('*.fhr');
curD=cd;
if length(d)>length(curD) && strcmp(d(1:length(curD)),curD)
    d=d(length(curD)+2:end);
end
set(handles.edit1,'String',[d,f]);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:length( handles.listPlotter )
    try delete( handles.listPlotter(i) ); catch, end
end
handles.listPlotter=fhrmorphoopenandplot( get(handles.edit1,'String'),...
    handles.listmethod(cell2mat(get(handles.uiMethodToggle(:),'Value'))==1,:) );
guidata(hObject, handles);
    
