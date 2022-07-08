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


classdef fhrplotFS < fhrplot
 

    properties
        EdtRec
        BtnNextRec
        BtnPrevRec
    end

    methods
        function obj=fhrplotFS(Sigs,TOCO,title,Sels)
            obj=obj@fhrplot(Sigs,TOCO,title,Sels);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function addSpecificControls(obj)
            obj.BtnDispPeriods=uicontrol(obj.Fig,'Style','toggle','String','Periodes','units','Pixels','position',[720 10 60 30],'Value',1,'Callback',@(src,evt) redraw(obj),'ForegroundColor',[.2 .2 .2]);            
            obj.BtnTgSigs(1)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{1},'units','Pixels','position',[780 10 80 30],'Value',1,'Callback',@(src,evt) redraw(obj));

            uicontrol(obj.Fig,'Style','text','String','MHR:','units','Pixels','position',[940 22 50 25]);
            uicontrol(obj.Fig,'Style','text','String','FHR:','units','Pixels','position',[1140 22 50 25]);
            for i=2:length(obj.SigNames)
                %FHR2,MHRe,MHRvale,MHR, data(2,:),FHR1e,FHR1vale,FHR1,data(1,:)
                val=(i==3 || i==5 || i==7 ||i==9 );
                if i<=5
                    obj.BtnTgSigs(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{i},'units','Pixels','position',[760+i*50 10 50 20],'Value',val);
                else
                    obj.BtnTgSigs(i)=uicontrol(obj.Fig,'Style','toggle','String',obj.SigNames{i},'units','Pixels','position',[765+i*50 10 50 20],'Value',val);
                end
                if (i==5 || i==9)
                set(obj.BtnTgSigs(i),'Callback',@(src,evt) redraw(obj))
                elseif i<=4
                    set(obj.BtnTgSigs(i),'Callback',@(src,evt) changeMHRDisp(obj,src))
                else
                    set(obj.BtnTgSigs(i),'Callback',@(src,evt) changeFHRDisp(obj,src))
                end
            end

            obj.EdtRec=uicontrol(obj.Fig,'Style','edit','String','','units','Pixels','position',[1310 10 50 30],'BackgroundColor',[1 1 1],'Callback',@(src,evt) changeRec(obj,str2double(get(obj.EdtRec,'String')))  );
            obj.BtnPrevRec=uicontrol(obj.Fig,'Style','pushbutton','String','< (P)','units','Pixels','position',[1270 10 40 30],'Callback',@(src,evt) changeRec(obj,str2double(get(obj.EdtRec,'String'))-1) );
            obj.BtnNextRec=uicontrol(obj.Fig,'Style','pushbutton','String','> (N)','units','Pixels','position',[1360 10 40 30],'Callback',@(src,evt) changeRec(obj,str2double(get(obj.EdtRec,'String'))+1) );
            obj.EdtComment=uicontrol(obj.Fig,'Style','edit','String','','units','Pixels','position',[1410 10 500 30],'BackgroundColor',[1 1 1],'Callback',@() obj.setNeedSave());

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function setConstants(obj) 
            obj.SelectTypeFile = 'SelectTypesFS.mat';
            obj.HelpImage='help_FS.png';
            obj.SigNames={'2nd Sensor','Raw','Gradient','Thresh 0.5','Interp','Raw','Gradient','Thresh 0.5','Interp'};
            obj.FHRRange=[30 255];
            obj.SigColors={[0 0 1;1 0 0;.6 .6 1;1 .6 .6],...
                       [.1 .6 .1],...
                       [1 0 1;0 .6 .6;1 .6 1;.3 .85 .85]};            
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeMHRDisp(obj,src)
            for i=2:4
                if obj.BtnTgSigs(i)~=src
                    set(obj.BtnTgSigs(i),'value',0)
                end
            end
            redraw(obj)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function changeFHRDisp(obj,src)
            for i=6:8
                if obj.BtnTgSigs(i)~=src
                    set(obj.BtnTgSigs(i),'value',0)
                end
            end
            redraw(obj)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function drawSigs(obj)
            d=floor(obj.Time*240)+1;
            f=min(size(obj.Sigs{1},2),floor(obj.Time*240+obj.WinLength*240));


            %FHR2
            if(get(obj.BtnTgSigs(1),'Value'))
                line(1/240:1/240:(f-d+1)/240,obj.Sigs{2}(1,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{2});
            end
            
            %Sig order:FHRraw; PFS; FHR.5; FHRrawInterp; FHR.5Interp

            for i=1:2 %[MHR FHR]
                chan=5-2*i; %[3 1]
                if(get(obj.BtnTgSigs(4*i+1),'Value')) %Interp (Btns 5 9)
                    if(get(obj.BtnTgSigs(4*i-2),'Value')) %Raw
                        line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(4,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(3,:)); % Display RawInter(4) in color Interp TS(3)
                    elseif(get(obj.BtnTgSigs(4*i-1),'Value') || get(obj.BtnTgSigs(4*i),'Value')) %Color gradient or Threshold at 0.5
                        line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(4,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(4,:)); % Display RawInter(4) in color Interp FS(4)
                        line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(5,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(3,:)); % Display ThreInter(5) in color Interp TS(3)
                    end
                end
                if(get(obj.BtnTgSigs(4*i-2),'Value')) %Raw (Btns 2 6)
                    line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(1,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(1,:)); %Raw in color TS
                elseif(get(obj.BtnTgSigs(4*i-1),'Value')) %Color gradient
                    P=obj.Sigs{chan}(2,d:f);
                    C=P'*obj.SigColors{chan}(2,:)+(1-P')*obj.SigColors{chan}(1,:);
                    x = 1/240:1/240:(f-d+1)/240;
                    y = obj.Sigs{chan}(1,d:f);
                    z = 3*ones(size(x));
                    col=reshape(C,[1 size(C)]);
                    surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','parent',obj.Axes);
                elseif(get(obj.BtnTgSigs(4*i),'Value')) %Theshold 0.5
                    line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(1,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(2,:)); %Raw(1) in color FS(2)
                    line(1/240:1/240:(f-d+1)/240,obj.Sigs{chan}(3,d:f),3*ones(1,f-d+1),'parent',obj.Axes,'Color',obj.SigColors{chan}(1,:)); %TS(3) in color TS(1)
                end
            end
            

            line(1/240:1/240:(f-d+1)/240,obj.TOCO(d:f),'parent',obj.AxesToco,'Color',[0 0 0]);

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseDown(obj)
            pos=get(obj.Axes,'CurrentPoint');
            ylim=get(obj.Axes,'Ylim');
            xlim=obj.WinLength;
            if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2)
                if (strcmp(obj.Mode,'PeriodSelect')) %Convert to rectangular selection
                    obj.SelectionStart=[pos(1,1)+obj.Time pos(1,2)];
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function MouseUp(obj)
            pos=get(obj.Axes,'CurrentPoint');
            ylim=get(obj.Axes,'Ylim');
            xlim=obj.WinLength;
            if (strcmp(obj.Mode,'PeriodSelect'))
                if pos(1,1)>=0 && pos(1,1)<=xlim && pos(1,2)>=ylim(1) && pos(1,2)<=ylim(2)
                    nkey=obj.KeySelection;
                    if nkey<=2, nSig=1; else, nSig=3;end
                    
                    D=min(pos(1,1)+obj.Time,obj.SelectionStart(1));
                    F=max(pos(1,1)+obj.Time,obj.SelectionStart(1));
                    S=obj.Sigs{nSig}(1,round(D*240+1):round(F*240));
                    V=S<max([pos(1,2) obj.SelectionStart(2)])& S>min([pos(1,2) obj.SelectionStart(2)]);

                    newSelect=(  round(D*240) + [ find(V&~[false V(1:end-1)])' find(V&~[V(2:end) false])' ]  )/240;
                    
                    if isempty(obj.Selections{nkey}), obj.Selections{nkey}=zeros(0,2); end
                    obj.Selections{nkey}=sortrows([obj.Selections{nkey};newSelect]);
                    i=1;
                    while i<size(obj.Selections{nkey},1)
                        interlaps=obj.Selections{nkey}(:,1)<=obj.Selections{nkey}(i,2) &...
                            obj.Selections{nkey}(:,2)>=obj.Selections{nkey}(i,1);
                        if sum(interlaps)>1
                            D=min(obj.Selections{nkey}(interlaps,1));
                            F=max(obj.Selections{nkey}(interlaps,2));
                            obj.Selections{nkey}=sortrows([obj.Selections{nkey}(~interlaps,:);[D F]]);
                            i=1;
                        else
                            i=i+1;
                        end
                    end
                    
                    obj.SelectionStart=[];
                    obj.KeySelection=[];
                    obj.NeedSave=true;
                end
                obj.Mode='';
                redraw(obj);
            end
        end        
    end
end