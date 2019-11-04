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


classdef fhrplot < hgsetget
    events
        Select
        ChangeTime
    end
    properties
        Fig
        Axes
        AxesToco
        BtnPrevious
        BtnPrevious1
        BtnNext1
        BtnNext
        BtnSelect
        BtnBaseLine
        BtnHist
        BtnSurf
        BtnPts
        EdtExp
        EdtInfo
        SelectionStart
        SelRectTmp
        LineCursor
        TextCursor
        chkbl
        BLnames
        BLcolors
        AccidentInfo
        BLexpert
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
        Selection
        SelectionAcc
        BaseLine
        BLPoints
        StatText
    end
    methods
        
        
        function changeTime(obj,src)
            if(src==obj.BtnPrevious)
                obj.Time=round(max(0,obj.Time-obj.WinLength+2));
            elseif(src==obj.BtnPrevious1)
                obj.Time=max(0,obj.Time-1);
            elseif(src==obj.BtnNext1)
                obj.Time=round(min(size(obj.Sigs{1},2)/240,obj.Time+1));
            elseif(src==obj.BtnNext)
                obj.Time=round(min(size(obj.Sigs{1},2)/240,obj.Time+obj.WinLength-2));
            end
            redraw(obj);
            notify(obj,'ChangeTime');
        end
        
        function obj=fhrplot(Sigs,TOCO,title,Acc,BLPoints)
            try
                g=load('defaultpos.mat');
            catch
                g.pos=[0, 100, 1000, 800];
            end
            accidents=load('accidents');
            obj.AccidentInfo=accidents.accidents;
            obj.SelectionAcc={};
            obj.SelectionAcc(1:size(obj.AccidentInfo,1))={zeros(0,2)};
            obj.BLPoints=zeros(2,0);
            if nargin<3, title=''; end
            obj.Fig=figure('Units','pixels','position',g.pos,'CloseRequestFcn',@(src,evts) close(obj),'WindowButtonMotionFcn',@(src,evt) MouseMovement(obj),'WindowButtonDownFcn',@(src,evt) MouseDown(obj),'ResizeFcn',@(src,evt) resize(obj),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt),'Renderer','zbuffer','Name',title,'numbertitle', 'off','toolbar','none');
           
            if nargin>=4, obj.SelectionAcc=Acc; end
            if nargin>=5, obj.BLPoints=BLPoints; end
            obj.Axes=axes('parent',obj.Fig,'Units','Pixels','position',[0 400 1594 500],'Xtick',[],'Ytick',[]);
            obj.AxesToco=axes('parent',obj.Fig,'Units','Pixels','position',[0 100 1594 250],'Xtick',[],'Ytick',[]);
            obj.BtnPrevious=uicontrol(obj.Fig,'Style','pushbutton','String','<<','units','Pixels','position',[10 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnPrevious1=uicontrol(obj.Fig,'Style','pushbutton','String','<','units','Pixels','position',[100 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnNext1=uicontrol(obj.Fig,'Style','pushbutton','String','>','units','Pixels','position',[190 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            obj.BtnNext=uicontrol(obj.Fig,'Style','pushbutton','String','>>','units','Pixels','position',[280 10 80 30],'Callback',@(src,evt) changeTime(obj,src),'KeyPressFcn', @(src,evt) KeyPressed(obj,evt));
            
            
            obj.EdtInfo=uicontrol(obj.Fig,'Style','edit','String','','units','Pixels','position',[380 10 180 30],'BackgroundColor',[1 1 1]);
            obj.BtnSelect=uicontrol(obj.Fig,'Style','pushbutton','String','save Select','units','Pixels','position',[570 10 80 30],'Callback',@(src,evt) clickSelect(obj));
            uicontrol(obj.Fig,'Style','pushbutton','String','Help','units','Pixels','position',[660 10 40 30],'Callback',@(src,evt) clickHelp(obj));
            uicontrol(obj.Fig,'Style','pushbutton','String','Eval','units','Pixels','position',[700 10 40 30],'Callback',@(src,evt) clickEval(obj));
            
            try obj.BLnames=blnames;   catch, for i=2:length(Sigs), obj.BLnames{i-1}=['Sig ' num2str(i-1)]; end,end
            try obj.BLcolors=blcolors; catch, obj.BLcolors=[0 0 0;.7 .7 .7;0 0 0;0 0 1;.7 .7 0;0 1 1;1 0 1;0 0 1;0 1 0;1 1 0;0 1 1;1 0 1];  end
            
            
            
            obj.Time=0;

            obj.StatText='';
            obj.Mode='Points';
            obj.Sigs=Sigs;
            obj.TOCO=TOCO;
            
            obj.Selection=zeros(0,2);
            
            obj.BaseLine=[];
            
            obj.SelectionStart=[];
            obj.KeySelection=[];
            obj.NeedSave=false;
            obj.BLexpert=zeros(1,length(obj.Sigs{1}));

            
            obj.chkbl=uicontrol(obj.Fig,'Style','toggle','String','Analysis','units','Pixels','position',[750 10 60 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[1 0 0]);
            
%             for i=2:length(Sigs)
%                 obj.chkbl(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.BLnames{i-1},'ForegroundColor',obj.BLcolors(i+1,:),'units','Pixels','position',[710+i*50 10 50 30],'Value',1);
%                 set(obj.chkbl(i),'Callback',@(src,evt) redraw(obj))
%             end              
            resize(obj);
            
          
        end
        
        function clickSelect(obj)
            obj.NeedSave=false;
            notify(obj,'Select');
        end
        
        function clickHelp(obj) %#ok<MANU>
            figure('position',[0 100 275 280])
            axes('position',[0 0 1 1])
            C=imread('help_legend.png');
            image(C)
        end
        
        function clickEval(obj)
            position=get(obj.Fig,'position');
            figure('position',[position(1:2) 275 400],'Toolbar','none','MenuBar','none','Name','Evaluation against expert','numbertitle', 'off')
            uicontrol('units','normalized','style','text','position',[0 0 1 1],'String',obj.StatText,'HorizontalAlignment','Left')
            
        end
        
        function clickBaseline(obj)
            obj.BaseLine=zeros(1,length(obj.Sigs{1}));
            P=ones(1,length(obj.Sigs{1}));
            
            for i=1:size(obj.Selection,1)
                t=round(obj.Selection(i,1)*240):round(obj.Selection(i,2)*240);
                P(t)=1;
                if length(t)<240
                    obj.BaseLine(t)=mean(obj.Sigs{1}(t));
                    t0=1;
                else
                    obj.BaseLine(t)=butterfilt(obj.Sigs{1}(t),240,0,3);
                    %obj.BaseLine(t)=obj.Sigs{1}(t);
                    t0=60;
                    
                end
                xf=mean(obj.BaseLine(t(1:t0)));
                %xf=obj.BaseLine(t(1));
                
                if i>1
                    tprev=round(obj.Selection(i-1,2)*240-t0+1);
                    xprev=mean(obj.BaseLine(tprev:round(obj.Selection(i-1,2)*240)));
                    v=linspace(xprev,xf,t(t0)-tprev+1);
                    c=[linspace(0,1,t0) ones(1,length(v)-2*t0) linspace(1,0,t0)];
                else
                    tprev=1;
                    xprev=xf;
                    v=linspace(xprev,xf,t(t0)-tprev+1);
                    c=[ ones(1,length(v)-t0) linspace(1,0,t0)];
                end
                
                obj.BaseLine(tprev:t(t0))=(1-c).*obj.BaseLine(tprev:t(t0))+c.*v;
            end
            obj.BaseLine(t(end):end)=obj.BaseLine(t(end));
            distancecoef=[linspace(0,1,200) 1   linspace(1,0,200)];
            obj.BaseLine=medgliss(obj.BaseLine,distancecoef.^16,P,24);
            %obj.BaseLine=butterfilt(obj.BaseLine,240,0,0.4);
            redraw(obj);
        end
        
        function redraw(obj)
            cla(obj.Axes);
            cla(obj.AxesToco);
            pos=get(obj.Axes,'Position');
            
            %[60-210]
            PixelsByCm=pos(4)/7.5;
            obj.WinLength=pos(3)/PixelsByCm;
            
            t=10*(ceil(obj.Time/10):floor((obj.Time+obj.WinLength)/10));
            
            set(obj.Axes,'Ylim',[60 210],'Xlim',[0 obj.WinLength],'Xtick',t-obj.Time,'XtickLabel',t);
            set(obj.AxesToco,'Ylim',[0 100],'Xlim',[0 obj.WinLength]);
            redrawSelect(obj);
            
            
            for i=1:7
                line([0 obj.WinLength],60+i*[20 20],'parent',obj.Axes,'Color',[0.6 1 0.6]);
            end
            line([0 obj.WinLength],[160 160],'parent',obj.Axes,'Color',[0.6 1 0.6],'LineWidth',2);
            line([0 obj.WinLength],[110 110],'parent',obj.Axes,'Color',[0.6 1 0.6],'LineWidth',2);
            for i=0:5
                line([0 obj.WinLength],i*[20 20],'parent',obj.AxesToco,'Color',[0.6 1 0.6]);
            end
            for i=1:obj.WinLength
                if(rem(obj.Time+i,10)==0)
                    e=2;
                else
                    e=1;
                end
                line([i i],[60 210],'parent',obj.Axes,'Color',[0.6 1 0.6],'LineWidth',e);
                line([i i],[0 100],'parent',obj.AxesToco,'Color',[0.6 1 0.6],'LineWidth',e);
            end
            %color=[0 0 0;0 0 1;1 0 0;0 1 0;1 1 0;0 1 1;1 0 1;0 0 1;1 0 0;0 1 0;1 1 0;0 1 1;1 0 1];
            
            d=floor(obj.Time*240)+1;
            
            f=min(size(obj.Sigs{1},2),floor(obj.Time*240+obj.WinLength*240));
            
            for i=1:length(obj.Sigs)
                if(size(obj.Sigs{i},1)==1)
                    if(i<=2 || get(obj.chkbl,'Value'))
                        line(1/240:1/240:(f-d+1)/240,obj.Sigs{i}(d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.BLcolors(i+1,:));
                    end
                else
                    if(i==1 || get(obj.chkbl,'Value'))
                        
                        patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.Sigs{i}(1,d:f) obj.Sigs{i}(2,f:-1:d)],1*ones(1,2*(f-d+1)),obj.BLcolors(i+2,:), 'FaceAlpha', 0.80,'parent',obj.Axes,'EdgeColor','none');
                        patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.Sigs{i}(1,d:f) obj.Sigs{i}(3,f:-1:d)],1*ones(1,2*(f-d+1)),obj.BLcolors(i+3,:), 'FaceAlpha', 0.80,'parent',obj.Axes,'EdgeColor','none');
                        patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.Sigs{i}(1,d:f) obj.Sigs{i}(4,f:-1:d)],1*ones(1,2*(f-d+1)),obj.BLcolors(i+4,:), 'FaceAlpha', 0.80,'parent',obj.Axes,'EdgeColor','none');
                        %patch([1/240:1/240:(f-d+1)/240 (f-d+1)/240:-1/240:1/240],[obj.Sigs{i}(1,d:f) obj.Sigs{i}(5,f:-1:d)],1*ones(1,2*(f-d+1)),obj.BLcolors(i+5,:), 'FaceAlpha', 0.80,'parent',obj.Axes,'EdgeColor','none');
                        
                        line(1/240:1/240:(f-d+1)/240,obj.Sigs{i}(1,d:f),2*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.BLcolors(i+1,:));
                    end
                end
            end
            line(1/240:1/240:(f-d+1)/240,obj.TOCO(d:f),'parent',obj.AxesToco,'Color',[0 0 0]);
            if(~isempty(obj.BaseLine) && get(obj.chkbl,'Value') )
                line(1/240:1/240:(f-d+1)/240,obj.BaseLine(d:f),'parent',obj.Axes,'Color',obj.BLcolors(1,:));
            end
            
            if(~isempty(obj.BLPoints) && get(obj.chkbl,'Value') )
                line(obj.BLPoints(1,:)-obj.Time,obj.BLPoints(2,:),5*ones(size(obj.BLPoints(1,:))),'parent',obj.Axes,'Color',[1 0.6 0.4],'Marker','s','LineStyle','none');
                xx=obj.Time:1/240:min(obj.Time+obj.WinLength,length(obj.Sigs{1})/240);
                x=[0,obj.BLPoints(1,:),length(obj.Sigs{1})/240];
                y=obj.BLPoints(2,[1 1:end end]);
                try
                    yy=linearinterpolation(x,y,xx);
                    line(xx-obj.Time,yy,5*ones(size(yy)),'parent',obj.Axes,'Color',[0.8 0.3 0.1])
                catch
                end
                %X=d+1:d+20*240;paperfig(obj.Sigs{1}(1,X),yy(1:length(X)),[obj.Sigs{13}(1,X);obj.Sigs{11}(1,X);obj.Sigs{12}(1,X);obj.Sigs{4}(1,X)],{'Mongelli','Lu','Wrobel','Pardey'})
                %ex1 : SCA050 min 41 decelration
                %ex2 : SCA001 min 18 pattern D
                %ex3 : Sca020 min 48 bl change
                %ex0 bed03-13285
                %
                %x=d+1:d+20*240;bls=vertcat(obj.Sigs{:});bls=bls([3:8 10:13],x);paperfig(obj.Sigs{1}(1,x),yy(1:length(x)),[min(bls);max(bls)],{'Min of 11 methods' 'Max of 11 methods'})
                
                
            end
            
        end
        
        function toggleHist(obj)
            if(get(obj.BtnHist,'Value'))
                obj.Mode='Hist';
                set(obj.BtnSurf,'Value',0);
            else
                obj.Mode='Points';
                set(obj.SelRectTmp,'FaceColor',[.8 1 .8]);
            end
            
        end
        
        function toggleSurf(obj)
            if(get(obj.BtnSurf,'Value'))
                obj.Mode='Surf';
                set(obj.BtnHist,'Value',0);
                set(obj.SelRectTmp,'FaceColor',[1 .8 .8]);
            else
                obj.Mode='Points';
                set(obj.SelRectTmp,'FaceColor',[.8 1 .8]);
            end
            
        end
        
        function MouseMovement(obj)
            
            ylim=get(obj.Axes,'Ylim');
            pos=get(obj.Axes,'CurrentPoint');
            if(~isempty(obj.SelectionStart) && (strcmp(obj.Mode,'Select') || strcmp(obj.Mode,'Surf') || ~isempty(obj.KeySelection)))
                set(obj.SelRectTmp,'Position',[min(obj.SelectionStart-obj.Time,pos(1,1)) ylim(1) abs(pos(1,1)-obj.SelectionStart+obj.Time)+0.0001 ylim(2)]);
            end
            set(obj.LineCursor(1),'XData',[pos(1,1) pos(1,1)],'YData',[60 210])
            set(obj.LineCursor(2),'XData',[0 obj.WinLength],'YData',[pos(1,2) pos(1,2)])
            set(obj.LineCursor(3),'XData',[pos(1,1) pos(1,1)],'YData',[0 100])
            set(obj.LineCursor(4),'XData',[pos(1,1)+0.25 pos(1,1)+0.25],'YData',[pos(1,2)-15 pos(1,2)+15])
            set(obj.LineCursor(5),'XData',[pos(1,1) pos(1,1)+0.25],'YData',[pos(1,2)+15 pos(1,2)+15])
            set(obj.LineCursor(6),'XData',[pos(1,1) pos(1,1)+0.25],'YData',[pos(1,2)-15 pos(1,2)-15])
            t=pos(1,1)+obj.Time;
            set(obj.EdtInfo,'String',[num2str(floor(t),'%02d') 'min' num2str(floor(rem(t*60,60)),'%02d') 's ' num2str(floor(10*pos(1,2))/10,'%05.1f') 'bpm']);
        end
        
        function KeyPressed(obj,evt)
            pos=get(obj.Axes,'CurrentPoint');
            ylim=get(obj.Axes,'Ylim');
            xlim=obj.WinLength;
            if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2) && get(obj.chkbl(1),'Value')
                nkey=find(strcmp(evt.Key,obj.AccidentInfo(:,1)));
                if ~isempty(nkey)
                    
                    if(strcmp(evt.Modifier,'shift'))
                        n=find(obj.SelectionAcc{nkey}(:,1)<=pos(1,1)+obj.Time & obj.SelectionAcc{nkey}(:,2)>=pos(1,1)+obj.Time);
                        if(~isempty(n))
                            obj.SelectionAcc{nkey}=obj.SelectionAcc{nkey}([1:n-1 n+1:end],:);
                        else
                            obj.SelectionStart=[];
                        end
                    else
                        if(isempty(obj.SelectionStart))
                            obj.SelectionStart=pos(1,1)+obj.Time;
                            obj.KeySelection=nkey;
                            
                            set(obj.SelRectTmp,'FaceColor',obj.AccidentInfo{nkey,3})
                        else
                            newSelect=[min(pos(1,1)+obj.Time,obj.SelectionStart) max(pos(1,1)+obj.Time,obj.SelectionStart)];
                            if isempty(obj.SelectionAcc{nkey}), obj.SelectionAcc{nkey}=zeros(0,2); end
                            E=obj.SelectionAcc{nkey}(:,1)<=newSelect(2) & obj.SelectionAcc{nkey}(:,2)>=newSelect(1);
                            newSelect(1)=min([obj.SelectionAcc{nkey}(E,1);newSelect(1)]);
                            newSelect(2)=max([obj.SelectionAcc{nkey}(E,2);newSelect(2)]);
                            obj.SelectionAcc{nkey}=[obj.SelectionAcc{nkey}(~E,:);newSelect];
                            obj.SelectionAcc{nkey}=sortrows(obj.SelectionAcc{nkey});
                            obj.SelectionStart=[];
                            obj.KeySelection=[];
                            obj.NeedSave=true;
                        end
                    end
                end
            end
            redraw(obj);
        end
        
        function MouseDown(obj)
            pos=get(obj.Axes,'CurrentPoint');
            ylim=get(obj.Axes,'Ylim');
            xlim=obj.WinLength;
            if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2)
                if (strcmp(obj.Mode,'Select'))
                    if(strcmp(get(obj.Fig,'SelectionType'),'alt'))
                        n=find(obj.Selection(:,1)<=pos(1,1)+obj.Time & obj.Selection(:,2)>=pos(1,1)+obj.Time);
                        if(~isempty(n))
                            obj.Selection=obj.Selection([1:n-1 n+1:end],:);
                        else
                            obj.SelectionStart=[];
                        end
                    else
                        if(isempty(obj.SelectionStart))
                            obj.SelectionStart=pos(1,1)+obj.Time;
                            
                        else
                            newSelect=[min(pos(1,1)+obj.Time,obj.SelectionStart) max(pos(1,1)+obj.Time,obj.SelectionStart)];
                            E=obj.Selection(:,1)<=newSelect(2) & obj.Selection(:,2)>=newSelect(1);
                            newSelect(1)=min([obj.Selection(E,1);newSelect(1)]);
                            newSelect(2)=max([obj.Selection(E,2);newSelect(2)]);
                            obj.Selection=[obj.Selection(~E,:);newSelect];
                            obj.Selection=sortrows(obj.Selection);
                            obj.SelectionStart=[];
                        end
                    end
                elseif(strcmp(obj.Mode,'Hist'))
                    expon=str2double(get(obj.EdtExp,'String'));
                    [histo,xi]=SamHisto(obj.Sigs{1},expon,round(240*(pos(1,1)+obj.Time)));
                    figure;
                    plot(xi,histo)
                elseif(strcmp(obj.Mode,'Surf'))
                    if(isempty(obj.SelectionStart))
                        obj.SelectionStart=pos(1,1)+obj.Time;
                        
                    else
                        srate=5;
                        expon=str2double(get(obj.EdtExp,'String'));
                        select=[min(pos(1,1)+obj.Time,obj.SelectionStart) max(pos(1,1)+obj.Time,obj.SelectionStart)];
                        
                        
                        T=select(1):1/srate:select(2);
                        SurfaceHist=zeros(length(T),180);
                        
                        for i=1:length(T)
                            [SurfaceHist(i,:),xi]=SamHisto(obj.Sigs{1},expon,round(select(1)*240+(i-1)*240/srate));
                        end
                        
                        
                        obj.SelectionStart=[];
                        figure;
                        surf(xi,T,SurfaceHist);
                    end
                elseif(strcmp(obj.Mode,'Points') && get(obj.chkbl(1),'Value'))
                    
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
%                     oldbl=obj.BLexpert;
%                     
%                     xx=1/240:1/240:length(obj.Sigs{1})/240;
%                     x=[0,obj.BLPoints(1,:),length(obj.Sigs{1})/240];
%                     y=obj.BLPoints(2,[1 1:end end]);
%                     obj.BLexpert=linearinterpolation(x,y,xx);
%                     t=(abs(oldbl-obj.BLexpert)>0.1);
%                     
%                     for i=1:length(obj.SelectionAcc)
%                         c=ones(1,size(obj.SelectionAcc{i},1));
%                         for j=1:size(obj.SelectionAcc{i},1)
%                             if t(round(240*obj.SelectionAcc{i}(j,1))) || t(round(240*obj.SelectionAcc{i}(j,2)))
%                                 c(j)=0;
%                             end
%                         end
%                         obj.SelectionAcc{i}=obj.SelectionAcc{i}(c==1,:);
%                     end
%                     
%                     [accelerations,decelerations]=accidentfromBL(obj.Sigs{1},obj.BLexpert);
%                     
%                     if ~isempty(accelerations)
%                         for j=1:size(accelerations,2)
%                             if t(round(4*accelerations(1,j))) || t(round(4*accelerations(2,j)))
%                                 obj.SelectionAcc{1}=[obj.SelectionAcc{1} ; accelerations(1:2,j)'/60];
%                             end
%                         end
%                     end
%                     if ~isempty(decelerations)
%                         for j=1:size(decelerations,2)
%                             if t(round(4*decelerations(1,j))) || t(round(4*decelerations(2,j)))
%                                 obj.SelectionAcc{2}=[obj.SelectionAcc{2} ; decelerations(1:2,j)'/60];
%                             end
%                         end
%                     end
%                     for j=1:2
%                         
%                         obj.SelectionAcc{j}=sortrows(obj.SelectionAcc{j},1);
%                         i=1;
%                         while i<size( obj.SelectionAcc{j},1)
%                             if obj.SelectionAcc{j}(i,2)>obj.SelectionAcc{j}(i+1,1)
%                                 obj.SelectionAcc{j}(i,2)=obj.SelectionAcc{j}(i+1,2);
%                                 obj.SelectionAcc{j}=obj.SelectionAcc{j}([1:i i+2:end],:);
%                             end
%                             i=i+1;
%                         end
%                     end
                    
                end
                
            end
            redraw(obj);
        end
        
        function redrawSelect(obj)
            ylim=get(obj.Axes,'Ylim');
            if(strcmp(obj.Mode,'Select'))
                cl=[.8 1 .8];
            else
                cl=[1 .8 .8];
            end
            if(get(obj.chkbl,'Value'))
                if(~isempty(obj.KeySelection))
                    cl=obj.AccidentInfo{obj.KeySelection,3};
                end
                obj.SelRectTmp=rectangle('Parent',obj.Axes,'Position',[-1 -1 0.001 0.001],'EdgeColor','none','FaceColor',cl);
                for i=1:size(obj.Selection,1)
                    if ( obj.Selection(i,1)<obj.Time+obj.WinLength && obj.Selection(i,2)>obj.Time)
                        rectangle('Parent',obj.Axes,'Position',[obj.Selection(i,1)-obj.Time ylim(1) obj.Selection(i,2)-obj.Selection(i,1)+0.00001 ylim(2)],'EdgeColor','none','FaceColor',[.8 1 .8]);
                    end
                end
                for j=1:size(obj.AccidentInfo,1)
                    for i=1:size(obj.SelectionAcc{j},1)
                        if ( obj.SelectionAcc{j}(i,1)<obj.Time+obj.WinLength && obj.SelectionAcc{j}(i,2)>obj.Time)
                            rectangle('Parent',obj.Axes,'Position',[obj.SelectionAcc{j}(i,1)-obj.Time ylim(1) obj.SelectionAcc{j}(i,2)-obj.SelectionAcc{j}(i,1)+0.00001 ylim(2)],'EdgeColor','none','FaceColor',obj.AccidentInfo{j,3});
                        end
                    end
                end
            end
            
            
            obj.LineCursor(1)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(2)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(4)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(5)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(6)=line('Parent',obj.Axes,'Color',[0 .5 0]);
            obj.LineCursor(3)=line('Parent',obj.AxesToco,'Color',[0 .5 0]);
            
            
        end
        
        function resize(obj)
            set(obj.Fig,'Units','pixels')
            pos=get(obj.Fig,'position');
            if(~isempty(pos))
                save('defaultpos.mat','pos');
                h=pos(4)-75;
                set(obj.Axes,'Position',[0 75+h/3 pos(3) h*2/3]);
                set(obj.AxesToco,'Position',[0 50 pos(3) h/3]);
                redraw(obj);
            end
        end
        
        function close(obj)
            if(obj.NeedSave)
                q=questdlg('Sauvegarder avant de quitter ?');
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