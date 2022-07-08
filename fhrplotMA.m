% Class which display an FHR viewer and enable expert to place baseline and
% to draw accidents
%
% USAGE 
%     obj=fhrplot(Sigs,TOCO)
%     obj=fhrplot(Sigs,TOCO,title)
%     obj=fhrplot(Sigs,TOCO,title,Acc)
%     obj=fhrplot(Sigs,TOCO,title,Acc,BLPoints)
%
% INPUT
%     Sigs , a Cell of row signal, for the top part generally {FHR,FHRraw,baseline}
%     TOCO the TOCOgraph sig
%     title  The figure title
%     Acc    cell list of periods {Acceleration, Deceleration, overShoots,Unreliable, Not to Analyse}
%     BLPoints List of expert baseline points
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


classdef fhrplotMA < fhrplot
    properties
        BtnDispXpertBaseLine
    end
    properties (Access=public)
        BLPoints
        StatText
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        function obj=fhrplotMA(Sigs,TOCO,title,Sels,BLPoints)
            obj=obj@fhrplot(Sigs,TOCO,title,Sels);
            if nargin<=4
                obj.BLPoints=zeros(2,0);
            else
                obj.BLPoints=BLPoints;
            end
            obj.redraw();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function setConstants(obj) 
            obj.SelectTypeFile = 'SelectTypesMA.mat';
            obj.HelpImage='help_MA.png';
            obj.SigNames={'Interpolation','MHR','Method'};
            obj.SigColors={[0 0 1],[.7 .7 1],[1 0 1]...
                          [1 0 0]}; 
            
        end  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function clickEval(obj)
            position=get(obj.Fig,'position');
            figure('position',[position(1:2) 275 400],'Toolbar','none','MenuBar','none','Name','Evaluation against expert','numbertitle', 'off')
            uicontrol('units','normalized','style','text','position',[0 0 1 1],'String',obj.StatText,'HorizontalAlignment','Left')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function addSpecificControls(obj)
            uicontrol(obj.Fig,'Style','pushbutton','String','Eval','units','Pixels','position',[700 10 40 30],'Callback',@(src,evt) clickEval(obj));
            obj.BtnDispPeriods=uicontrol(obj.Fig,'Style','toggle','String','Expert Per.','units','Pixels','position',[750 10 60 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[.2 .2 .2]);  
            obj.BtnDispXpertBaseLine=uicontrol(obj.Fig,'Style','toggle','String','Expert BL','units','Pixels','position',[815 10 60 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[.2 .2 .2]);  
            
            obj.BtnTgSigs(1)=uicontrol(obj.Fig,'Style','toggle','String','FHR','ForegroundColor',obj.SigColors{1},'units','Pixels','position',[-100 10 50 30],'Value',1,'Visible','off','Callback',@(src,evt) redraw(obj));
            for i=2:length(obj.Sigs)
                obj.BtnTgSigs(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{i-1},'ForegroundColor',obj.SigColors{i},'units','Pixels','position',[785+i*50 10 50 30],'Value',1,'Callback',@(src,evt) redraw(obj));
            end 
            for i=4:length(obj.Sigs)
                obj.BtnTgSigs(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{i-1},'ForegroundColor',obj.SigColors{i},'units','Pixels','position',[785+i*50 10 50 30],'Value',1,'Callback',@(src,evt) changeMethod(obj,src));
            end             
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeMethod(obj,src)
            for i=4:length(obj.BtnTgSigs)
                if obj.BtnTgSigs(i)~=src
                    set(obj.BtnTgSigs(i),'value',0)
                end
            end
            redraw(obj)            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function drawSigs(obj)
            obj.drawSigs@fhrplot();
            if(~isempty(obj.BLPoints) && get(obj.BtnDispXpertBaseLine,'Value') )
                line(obj.BLPoints(1,:)-obj.Time,obj.BLPoints(2,:),5*ones(size(obj.BLPoints(1,:))),'parent',obj.Axes,'Color',[1 0.6 0.4],'Marker','s','LineStyle','none');
                xx=obj.Time:1/240:min(obj.Time+obj.WinLength,length(obj.Sigs{1})/240);
                x=[0,obj.BLPoints(1,:),length(obj.Sigs{1})/240];
                y=obj.BLPoints(2,[1 1:end end]);
                %try
                    yy=linearinterpolation(x,y,xx);
                    line(xx-obj.Time,yy,5*ones(size(yy)),'parent',obj.Axes,'Color',[0.8 0.3 0.1])
                %catch
                %end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseDown(obj)
            pos=get(obj.Axes,'CurrentPoint');
            ylim=get(obj.Axes,'Ylim');
            xlim=obj.WinLength;
            if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2)
                if(strcmp(obj.Mode,'') && get(obj.BtnDispXpertBaseLine,'Value'))
                    if(strcmp(get(obj.Fig,'SelectionType'),'alt'))
                        n=find(obj.BLPoints(1,:)>=pos(1,1)+obj.Time-.2 & obj.BLPoints(1,:)<=pos(1,1)+obj.Time+.2 & ...
                            obj.BLPoints(2,:)>=pos(1,2)-2 & obj.BLPoints(2,:)<=pos(1,2)+2);
                        if(~isempty(n))
                            obj.BLPoints=obj.BLPoints(:,[1:n-1 n+1:end]);
                            obj.NeedSave=true;
                        end
                    else
                        if isempty(obj.BLPoints) || min(abs(pos(1,1)+obj.Time-obj.BLPoints(1,:))>0.1)
                            obj.BLPoints=[obj.BLPoints [pos(1,1)+obj.Time;pos(1,2)]];
                            obj.BLPoints=sortrows(obj.BLPoints')';
                            obj.NeedSave=true;
                        end
                    end
                end
            end
            redraw(obj);
        end     
    end
end