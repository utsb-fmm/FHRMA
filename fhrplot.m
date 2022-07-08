% Simple FHR viewer. Can Display FHR1,FHR2,MHR,TOCO and an FHR baseline
%
% fhrplotFS is inherited from specific tools for display and annotate False Signals 
% fhrplotMA is inherited from specific tools for display and annotate Morphological analysis (baseline, Acceleration, Deceleration)
%
% USAGE 
%     obj=fhrplot(FHRFile)
%     obj=fhrplot(Sigs,TOCO)
%     obj=fhrplot(Sigs,TOCO,title)
%     obj=fhrplot(Sigs,TOCO,title,Sels)
%
% INPUT
%     FHRFile a path to a .fhr, .fhrm or .dat file
%     Sigs , a Cell of row signal, expected {FHR1,FHR2,MHR,BaselineFHR} but can be only 1 signal
%     TOCO the TOCOgraph sig
%     title  The figure title
%     Sels a cell of n matrices size (ni x 2) indicating start and end of
%          each period of Selection on a row. There is n types of Selection each in different color.  
%
%
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


classdef fhrplot < matlab.mixin.SetGet
    events
        Save
        ChangeTime
        ChangeRec
    end
    properties
        Fig
        Axes
        AxesToco
        BtnPrevious
        BtnPrevious1
        BtnNext1
        BtnNext
        BtnSave
        EdtInfo
        EdtComment
        BtnDispPeriods
        BtnTgSigs
        BtnCU
        SelectionStart
        SelRectTmp
        LineCursor
        TextCursor
        SelectInfos
    end
    properties (SetAccess=protected) %Readonly properties
        Sigs
        TOCO
        Mode
        NeedSave
        KeySelection
    end
    properties (SetAccess=public)
        Time
        WinLength
        Selections
    end
    properties (Access = protected)
        SelectTypeFile = 'SelectTypesStd.mat'
        HelpImage='help_std.png'
        SigNames={'FHR1','FHR2','MHR','Baseline&A/D','FHR interp','RawFHR1','RawFHR2','RawMHR'}
        SigColors={[0 0 1],[0 0 1],[1 .4 1],[.1 .6 .1],[.7 .7 .7],[.9 .65 .65],[.5 .5 1],[1 .5 1]}
        PatchColors=[1 .5 .5;.5 1 .5;.5 .5 1];        
        FHRRange=[50 210];
        cmPerMinutes=1;
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj=fhrplot(Sigs,TOCO,title,Sels)
            if nargin==1
                [~,fname,ext]=fileparts(Sigs);
                title=[fname,ext];
                [FHR1,FHR2,MHR,TOCO]=fhropen(Sigs);
                Sigs={FHR1,FHR2,MHR};
            elseif nargin<3
                title=''; 
            end
            
            obj.setConstants()

            try
                g=load('defaultpos.mat');
            catch
                g.pos=[0, 100, 1000, 800];
            end
            
            obj.Time=0;
            obj.Mode='';
            obj.Sigs=Sigs;
            for i=1:length(obj.Sigs)
                obj.Sigs{i}(obj.Sigs{i}(:)==0)=NaN;
            end
            obj.TOCO=TOCO;

            obj.SelectionStart=[];
            obj.KeySelection=[];
            obj.NeedSave=false;

            
            obj.Fig=figure('Units','pixels','position',g.pos,'CloseRequestFcn',@(src,evts) close(obj),...
                'WindowButtonMotionFcn',@(src,evt) MouseMovement(obj),...
                'WindowButtonDownFcn',@(src,evt) MouseDown(obj),...
                'WindowButtonUpFcn',@(src,evt) MouseUp(obj),...
                'ResizeFcn',@(src,evt) resize(obj),...
                'KeyPressFcn', @(src,evt) KeyPressed(obj,evt),...
                'Renderer','zbuffer','Name',title,'numbertitle', 'off','toolbar','none','menubar','none');

            SelectTypes=load(obj.SelectTypeFile);
            obj.SelectInfos=SelectTypes.SelectTypes;
            obj.Selections={};
            obj.Selections(1:size(obj.SelectInfos,1))={zeros(0,2)};
            if nargin>=4, obj.Selections=Sels; end

            obj.Axes=axes('parent',obj.Fig,'Units','Pixels','position',[0 400 1594 500],'Xtick',[],'Ytick',[],'FontSize',14);
            obj.AxesToco=axes('parent',obj.Fig,'Units','Pixels','position',[0 100 1594 250],'Xtick',[],'Ytick',[]);
            obj.BtnPrevious=uicontrol(obj.Fig,'Style','pushbutton','String','<< (X)','units','Pixels','position',[10 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnPrevious1=uicontrol(obj.Fig,'Style','pushbutton','String','< (C)','units','Pixels','position',[100 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnNext1=uicontrol(obj.Fig,'Style','pushbutton','String','> (V)','units','Pixels','position',[190 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnNext=uicontrol(obj.Fig,'Style','pushbutton','String','>> (B)','units','Pixels','position',[280 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            
            obj.EdtInfo=uicontrol(obj.Fig,'Style','edit','String','','units','Pixels','position',[380 10 180 30],'BackgroundColor',[1 1 1]);
            obj.BtnSave=uicontrol(obj.Fig,'Style','pushbutton','String','Save','units','Pixels','position',[570 10 80 30],'Callback',@(src,evt) clickSave(obj));
            uicontrol(obj.Fig,'Style','pushbutton','String','Help','units','Pixels','position',[660 10 40 30],'Callback',@(src,evt) clickHelp(obj));
            obj.addSpecificControls();           
 
            resize(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function setConstants(obj)        %#ok<MANU> 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function addSpecificControls(obj)
            obj.BtnDispPeriods=uicontrol(obj.Fig,'Style','toggle','String','Periodes','units','Pixels','position',[720 10 60 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[.2 .2 .2]);
            if size(obj.TOCO,1)>1
                obj.BtnCU=uicontrol(obj.Fig,'Style','toggle','String','UC','units','Pixels','position',[785 10 40 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[0 .5 0]);
            end
            for i=1:length(obj.Sigs)
                obj.BtnTgSigs(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{i},'ForegroundColor',obj.SigColors{i},'units','Pixels','position',[750+i*80 10 80 30],'Value',1,'Callback',@(src,evt) redraw(obj));
            end 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeTime(obj,src)
            if(src==obj.BtnPrevious)
                obj.Time=round(max(0,obj.Time-obj.WinLength+2));
            elseif(src==obj.BtnPrevious1)
                obj.Time=max(0,obj.Time-1);
            elseif(src==obj.BtnNext1)
                obj.Time=round(min(size(obj.Sigs{1},2)/240-1,obj.Time+1));
            elseif(src==obj.BtnNext)
                obj.Time=round(min(size(obj.Sigs{1},2)/240-1,obj.Time+obj.WinLength-2));
            end
            redraw(obj);
            notify(obj,'ChangeTime');
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setNeedSave(obj)
            obj.NeedSave=true;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clickSave(obj)
            obj.NeedSave=false;
            notify(obj,'Save');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clickHelp(obj) 
            C=imread(obj.HelpImage);
            figure('position',[0 100 size(C,2) size(C,1)])
            axes('position',[0 0 1 1])
            image(C)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function redraw(obj)
            %*********Axes*********
            cla(obj.Axes);
            cla(obj.AxesToco);
            pos=get(obj.Axes,'Position');
            
            PixelsBySquare=20*pos(4)/(obj.FHRRange(2)-obj.FHRRange(1));
            obj.WinLength=pos(3)/PixelsBySquare/obj.cmPerMinutes;
            
            t=5*(ceil(obj.Time/5):floor((obj.Time+obj.WinLength)/5));
            
            set(obj.Axes,'Ylim',obj.FHRRange,'Xlim',[0 obj.WinLength],'Xtick',t-obj.Time,'XtickLabel',t);
            set(obj.AxesToco,'Ylim',[0 100],'Xlim',[0 obj.WinLength]);
            
            rectangle('Position',[0 110 obj.WinLength 50],'Parent',obj.Axes,'FaceColor',[.97 .97 .97],'LineStyle','none')

            %*********Periods*********
            redrawSelect(obj);
            
            %*********Grid*********

            T10=find(rem(obj.Time+(1:obj.WinLength),10)==0);
            squareWidth=18/PixelsBySquare/obj.cmPerMinutes;
            Xlines=[T10-squareWidth;NaN*ones(size(T10));T10+squareWidth];
            Xlines=[0 Xlines(:)' obj.WinLength];
            Ylines=[ (60:20:200)-squareWidth*15; NaN*ones(1,8); (60:20:200)+squareWidth*15];
            Ylines=[50 Ylines(:)' 210];
            YlinesTOCO=[ (20:20:80)-squareWidth*15; NaN*ones(1,4); (20:20:80)+squareWidth*15];
            YlinesTOCO=[0 YlinesTOCO(:)' 100];            

            for i=.5:.25:8.25
                e=1+(rem(i,1)==0);
                a=1-2.5*rem(i,.5);
                if rem(i,1)==0
                    line(Xlines,40+i*20*ones(size(Xlines)),'parent',obj.Axes,'Color',[0.66 1 0.66 a],'LineWidth',e);
                else
                    line([0 obj.WinLength],40+i*20*[1 1],'parent',obj.Axes,'Color',[0.66 1 0.66 a],'LineWidth',e);
                end
            end
            line(Xlines,160*ones(size(Xlines)),'parent',obj.Axes,'Color',[0.66 1 0.66],'LineWidth',4);
            line([0 obj.WinLength],[110 110],'parent',obj.Axes,'Color',[0.66 1 0.66],'LineWidth',4);
            for i=0:.5:5
                e=1+2*(i==1);
                a=1-rem(i,1);
                if rem(i,1)==0
                    line(Xlines,i*20*ones(size(Xlines)),'parent',obj.AxesToco,'Color',[0.66 1 0.66 a],'LineWidth',e);
                else
                    line([0 obj.WinLength],i*[20 20],'parent',obj.AxesToco,'Color',[0.66 1 0.66 a],'LineWidth',e);
                end
            end
            for i=.5/obj.cmPerMinutes:.5/obj.cmPerMinutes:obj.WinLength

                e=2+(rem(obj.Time+i,5)==0)+(rem(obj.Time+i,10)==0);
                a=1-.6*(rem(obj.Time+i,1)~=0);
                if rem(obj.Time+i,10)==0
                    line(i*ones(size(Ylines)),Ylines,'parent',obj.Axes,'Color',[0.66 1 0.66 a],'LineWidth',e);
                    line(i*ones(size(YlinesTOCO)),YlinesTOCO,'parent',obj.AxesToco,'Color',[0.66 1 0.66 a],'LineWidth',e);
                else
                    line([i i],[50 210],'parent',obj.Axes,'Color',[0.66 1 0.66 a],'LineWidth',e);
                    line([i i],[0 100],'parent',obj.AxesToco,'Color',[0.66 1 0.66 a],'LineWidth',e);

                end
                if(rem(obj.Time+i,10)==0)
                    for j=0:7
                        text(i,60+j*20,0,num2str(60+j*20),...
                            'Parent',obj.Axes,...
                            'HorizontalAlignment','center',...
                            'Color',[0.66 1 0.66],...
                            'FontSize',14)
                    end
                    for j=1:4
                        text(i,j*20,0,num2str(j*20),...
                            'Parent',obj.AxesToco,...
                            'HorizontalAlignment','center',...
                            'Color',[0.66 1 0.66],...
                            'FontSize',14)
                    end
                end
            end
            %*********Heart Rate Sigs*********
            obj.drawSigs()
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function drawSigs(obj)
            d=floor(obj.Time*240)+1;
            f=min(size(obj.Sigs{1},2),floor(obj.Time*240+obj.WinLength*240));

            for i=length(obj.Sigs):-1:1
                if(get(obj.BtnTgSigs(i),'Value'))
                    for j=2:size(obj.Sigs{i},1)
                        patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.Sigs{i}(1,d:f) obj.Sigs{i}(j,f:-1:d)],1*ones(1,2*(f-d+1)),obj.PatchColors(j-1,:), 'FaceAlpha', 0.80,'parent',obj.Axes,'EdgeColor','none');
                    end
                    line(1/240:1/240:(f-d+1)/240,obj.Sigs{i}(1,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{i});
                end
            end
            
            if size(obj.TOCO,1)>1 && get(obj.BtnCU,'Value')
                patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.TOCO(1,d:f) obj.TOCO(2,f:-1:d)],1*ones(1,2*(f-d+1)),[.5 1 .2], 'FaceAlpha', 0.80,'parent',obj.AxesToco,'EdgeColor','none');
            end
            line(1/240:1/240:(f-d+1)/240,obj.TOCO(1,d:f),ones(1,f-d+1),'parent',obj.AxesToco,'Color',[0 0 0]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseMovement(obj)
            
            ylim=get(obj.Axes,'Ylim');
            pos=get(obj.Axes,'CurrentPoint');
            posTOCO=get(obj.AxesToco,'CurrentPoint');

            if(~isempty(obj.SelectionStart) && ~isempty(obj.KeySelection))
                if(length(obj.SelectionStart)==2)
                    set(obj.SelRectTmp,'Position',[min(obj.SelectionStart(1)-obj.Time,pos(1,1)) min(obj.SelectionStart(2),pos(1,2)) abs(pos(1,1)-obj.SelectionStart(1)+obj.Time)+0.0001 abs(obj.SelectionStart(2)-pos(1,2)) ]);
                else
                    set(obj.SelRectTmp,'Position',[min(obj.SelectionStart-obj.Time,pos(1,1)) ylim(1) abs(pos(1,1)-obj.SelectionStart+obj.Time)+0.0001 ylim(2)]);
                end
            end
            set(obj.LineCursor(1),'XData',[pos(1,1) pos(1,1)],'YData',ylim ,'Zdata',[1 1])
            set(obj.LineCursor(2),'XData',[0 obj.WinLength],'YData',[pos(1,2) pos(1,2)],'Zdata',[1 1])
            set(obj.LineCursor(3),'XData',[pos(1,1) pos(1,1)],'YData',[0 100],'Zdata',[1 1])
            set(obj.LineCursor(4),'XData',[pos(1,1)+0.25 pos(1,1)+0.25],'YData',[pos(1,2)-15 pos(1,2)+15],'Zdata',[1 1])
            set(obj.LineCursor(5),'XData',[pos(1,1) pos(1,1)+0.25],'YData',[pos(1,2)+15 pos(1,2)+15],'Zdata',[1 1])
            set(obj.LineCursor(6),'XData',[pos(1,1) pos(1,1)+0.25],'YData',[pos(1,2)-15 pos(1,2)-15],'Zdata',[1 1])
            set(obj.LineCursor(7),'XData',[0 obj.WinLength],'YData',[posTOCO(1,2) posTOCO(1,2)],'Zdata',[1 1])
            if(pos(1,2)>=ylim(1)-5)
                ytext=[num2str(floor(10*pos(1,2))/10,'%05.1f') ' bpm'];
            else
                ytext=[num2str(floor(10*posTOCO(1,2))/10,'%04.1f') ' mmHg'];
            end
            t=pos(1,1)+obj.Time;
            set(obj.EdtInfo,'String',[num2str(floor(t),'%02d') ' min ' num2str(floor(rem(t*60,60)),'%02d') ' s - ' ytext]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function KeyPressed(obj,evt)
            if length(evt.Key)==1 && (evt.Key=='n' || evt.Key=='N')
                changeRec(obj,str2double(get(obj.EdtRec,'String'))+1)
            elseif length(evt.Key)==1 && ( evt.Key=='p'|| evt.Key=='P')
                changeRec(obj,str2double(get(obj.EdtRec,'String'))-1)
            elseif length(evt.Key)==1 && ( evt.Key=='v' || evt.Key=='V')
                changeTime(obj,obj.BtnNext1)                
            elseif length(evt.Key)==1 && ( evt.Key=='b' || evt.Key=='B')
                changeTime(obj,obj.BtnNext)
            elseif length(evt.Key)==1 && ( evt.Key=='c' || evt.Key=='C')
                changeTime(obj,obj.BtnPrevious1)
            elseif length(evt.Key)==1 && ( evt.Key=='x' || evt.Key=='X')
                changeTime(obj,obj.BtnPrevious)                
            else

                pos=get(obj.Axes,'CurrentPoint');
                ylim=get(obj.Axes,'Ylim');
                xlim=obj.WinLength;
                if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2) && get(obj.BtnDispPeriods,'Value')
                    nkey=find(strcmp(evt.Key,obj.SelectInfos(:,1)));
                    if ~isempty(nkey) 
                        if(strcmp(evt.Modifier,'shift')) %Shift + letter remove  the Periods
                            n=find(obj.Selections{nkey}(:,1)<=pos(1,1)+obj.Time & obj.Selections{nkey}(:,2)>=pos(1,1)+obj.Time);
                            if(~isempty(n))
                                obj.Selections{nkey}=obj.Selections{nkey}([1:n-1 n+1:end],:);
                            else
                                obj.SelectionStart=[];
                            end
                        else
                            if(isempty(obj.SelectionStart)) % Press first time so start Selection
                                obj.SelectionStart=pos(1,1)+obj.Time;
                                obj.KeySelection=nkey;
                                set(obj.SelRectTmp,'FaceColor',obj.SelectInfos{nkey,3})
                                obj.Mode='PeriodSelect';
                            else % Press Second time so end selection
                                if length(obj.SelectionStart)==1 %No rectangular selection in progres
                                    newSelect=[min(pos(1,1)+obj.Time,obj.SelectionStart) max(pos(1,1)+obj.Time,obj.SelectionStart)];
                                    if isempty(obj.Selections{nkey}), obj.Selections{nkey}=zeros(0,2); end
                                    E=obj.Selections{nkey}(:,1)<=newSelect(2) & obj.Selections{nkey}(:,2)>=newSelect(1);
                                    newSelect(1)=min([obj.Selections{nkey}(E,1);newSelect(1)]);
                                    newSelect(2)=max([obj.Selections{nkey}(E,2);newSelect(2)]);
                                    obj.Selections{nkey}=[obj.Selections{nkey}(~E,:);newSelect];
                                    obj.Selections{nkey}=sortrows(obj.Selections{nkey});
                                end
                                obj.SelectionStart=[];
                                obj.KeySelection=[];
                                obj.NeedSave=true;
                                obj.Mode='';
                            end
                        end
                    end
                end
                redraw(obj);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseUp(obj) %#ok<MANU> 
            %To overide
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseDown(obj) %#ok<MANU> 
            %To overide
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function redrawSelect(obj)
            ylim=get(obj.Axes,'Ylim');

            if(get(obj.BtnDispPeriods,'Value'))
                if isempty(obj.KeySelection)
                    cl=[.8 .8 .8];
                else
                    cl=obj.SelectInfos{obj.KeySelection,3};
                end
                
                for j=size(obj.SelectInfos,1):-1:1
                    for i=1:size(obj.Selections{j},1)
                        if ( obj.Selections{j}(i,1)<obj.Time+obj.WinLength && obj.Selections{j}(i,2)>obj.Time)
                            rectangle('Parent',obj.Axes,'Position',[obj.Selections{j}(i,1)-obj.Time ylim(1)+(5-j)*2 obj.Selections{j}(i,2)-obj.Selections{j}(i,1)+0.00001 ylim(2)-ylim(1)-2*(5-j)*2],'EdgeColor','none','FaceColor',obj.SelectInfos{j,3});
                        end
                    end
                end
                obj.SelRectTmp=rectangle('Parent',obj.Axes,'Position',[-1 -1 0.001 0.001],'EdgeColor','none','FaceColor',cl);
            end
            
            %********************************
            obj.LineCursor(1)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(2)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(4)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(5)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(6)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(3)=line('Parent',obj.AxesToco,'Color',[0 .5 0]);  
            obj.LineCursor(7)=line('Parent',obj.AxesToco,'Color',[0 .5 0]);  

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function resize(obj)
            set(obj.Fig,'Units','pixels')
            pos=get(obj.Fig,'position');
            if(~isempty(pos))
                save('defaultpos.mat','pos');
                h=pos(4)-82;
                set(obj.Axes,'Position',[0 82+h/3 pos(3) h*2/3]);
                set(obj.AxesToco,'Position',[0 50 pos(3) h/3]);
                redraw(obj);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeRec(obj,val)
           set(obj.EdtRec,'String',val)
           notify(obj,'ChangeRec') 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function close(obj)
            if(obj.NeedSave)
                q=questdlg('Save before exit ?');
                if(strcmp(q,'Yes'))
                    notify(obj,'Select');
                end
                if(strcmp(q,'No') ||strcmp(q,'Yes'))
                    pause(0.2)
                    obj.delete();
                end
            else
                obj.delete();
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function delete(obj)
            h = obj.Fig;
            if ishandle(h)
                delete(h);
            else
                return
            end
        end
        
    end
end