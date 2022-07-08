% This function open an FHR file, compute FSDop or FSScalp and launch the
% viewer
%
% USAGE
%     plotter=showFSAnalysis(file,chan)
% INPUT 
%     file : the file path of the FHR file
%     chan : 'Dop' or 'Scalp' depending on the channel to analyse and the
%     model to apply for FS prediction.


%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     FHR Morphological Analysis Toolbox  Copyright (C) 2022 Samuel Boudet, Faculté de Médecine et Maïeutique,
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

function plotter=showFSAnalysis(file,chan,Stage2Start)

if nargin==0
    [file,folder]=uigetfile({'*.fhrm';'*.dat'},...
               'Select an FHR file','FSdataset');
    file=[folder file];
end
if nargin<=1
    chan=questdlg('Which FHR channel have to be analyzed?', ...
        'Choose FHR channel', ...
        'Doppler (FHR1)','Scalp ECG (FHR2)','Doppler (FHR1)');
    if strcmp(chan,'Doppler (FHR1)')
        chan='Dop';
    elseif strcmp(chan,'Scalp ECG (FHR2)')
        chan='Scalp';
    else
        exit();
    end
end

[FHR1,FHR2,MHR,TOCO]=fhropen(file);
isStage2=zeros(size(FHR1));
load('FSdataset\expertAnnotations.mat','DB');

[~,fname,ext]=fileparts(file);
fname=[fname,ext];
n=find(strcmp(fname,{DB.filename}));

if nargin<=3
    if isempty(n)
        s2=inputdlg(sprintf(['When does the second stage start (in min) ?\n' ...
            'Leave empty if no second stage.\n' ...
            'Use 0 if the recording is only second stage.\n' ...
            'For CTG-UHB database, the second stage start at 60 min for all recs.']));
        Stage2Start=str2double(s2{1})*240+1;
    else
        Stage2Start=DB(n).Stage2_Start;
    end
end
if ~isempty(Stage2Start) && ~isnan(Stage2Start) && isnumeric(Stage2Start)
    isStage2(Stage2Start:end)=1;
end

if strcmp(chan,'Dop') 
    [PFHR,PMHR,FHRtrue,MHRtrue]=FalseSigDetectDopMHR(FHR1,MHR,isStage2);
elseif strcmp(chan,'Scalp') 
    [PFHR,FHRtrue,FHR2]=FalseSigDetectScalp(FHR2,isStage2);
    PMHR=zeros(size(FHR1));MHRtrue=MHR;
    FHRtmp=FHR2;
    FHR2=FHR1;
    FHR1=FHRtmp;
else
    error('chan must be Dop or Scalp');    
end
FHR1(FHR1==0)=NaN;FHR2(FHR2==0)=NaN;MHR(MHR==0)=NaN;
FHRtrue(FHRtrue==0)=NaN;MHRtrue(MHRtrue==0)=NaN;

FHR1=[FHR1;PFHR;FHRtrue;interpolFHR(FHR1);interpolFHR(FHRtrue)];
MHR=[MHR;PMHR;MHRtrue;interpolFHR(MHR);interpolFHR(MHRtrue)];



if ~isempty(n)
    Selections={DB(n).TS_FHR/240, DB(n).FS_FHR/240, DB(n).TS_MHR/240, DB(n).FS_MHR/240, zeros(0,2)};
else
    Selections(1:5)={zeros(0,2)};
    n=length(DB)+1;
    DB(n)=struct('filename',fname,...
        'dataset',[chan 'Custom'],...
        'TS_FHR',zeros(0,2),...
        'FS_FHR',zeros(0,2),...
        'TS_MHR',zeros(0,2),...
        'FS_MHR',zeros(0,2),...
        'Stage2_Start',(Stage2Start-1)/240, ...
        'Comment','');
    DB(n).filename=fname;
    DB(n).dataset=[chan 'Custom'];
    %save('FSdataset\expertAnnotations.mat','DB');
end

plotter=fhrplotFS({FHR1,FHR2,MHR},TOCO, fname,Selections); 
set(plotter.EdtComment,'String',DB(n).Comment)
addlistener(plotter,'Save',@(src,evtdat) Save(fname,plotter.Selections,get(plotter.EdtComment,'String')));
end

function Save(fname,Selections,comment)
    load('FSdataset\expertAnnotations.mat','DB');
    n=find(strcmp(fname,{DB.filename}));
    DB(n).Comment=comment;
    DB(n).TS_FHR=Selections{1};
    DB(n).FS_FHR=Selections{2};
    DB(n).TS_MHR=Selections{3};
    DB(n).FS_MHR=Selections{4};
    save('FSdataset\expertAnnotations.mat','DB');
end