% Open an FHR file, compute the sevelal analyses and plot them
%
%
% USAGE
%    plotters=fhrmorphoopenandplot(filename,analyses)
% INPUT
%    filename           the fhr file path
%    analyses           cell array of the command to get the analysis
%                           ex 'aamlu(FHR)' ; 
%                             FHR is the interpollated FHR;
%                             FHR0 is the raw FHR with 0 on missing signal
%                             TOCO is the TOC
%                             Each command should return
%                             [baseline,acceleration,deceleration]
               

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



function plotters=fhrmorphoopenandplot(filename,analyses)
    [FHR1,FHR2,TOCO]=fhropen(filename);
    [FHR,FHRraw,TOCO]=preprocess(FHR1,FHR2,TOCO);
    FHR0=FHRraw;FHR0(isnan(FHR0))=0;
    g=load('expertAnalyses');
    ff=find(['\' filename]=='\'| ['\' filename]=='/',1,'last');
    s=find(strcmpi(filename(ff:end), {g.data(:).filename} ));
    if ~isempty(s)
        SelectionAcc={g.data(s).accelerations, g.data(s).decelerations, g.data(s).overshoots, g.data(s).unreliableSignal, g.data(s).notToAnalyse};
    else
        SelectionAcc(1:5)={zeros(0,2)};
    end
    for i=1:size(analyses,1)
        if strcmp(analyses{i,2},'expert') 
            if ~isempty(s), BLPoints=g.data(s).expertPts;else BLPoints=zeros(2,0);end
            plotters(i)=fhrplot({FHR,FHRraw},TOCO,[analyses{i,1} ' - ' filename],SelectionAcc,BLPoints); %#ok<*AGROW>
            
            
            addlistener(plotters(i),'Select',@(src,evtdat) Select(filename,plotters(i).BLPoints,plotters(i).SelectionAcc));
        else
            S=SelectionAcc;S{3}=zeros(0,2);
            [baseline,acc,dec]=eval(analyses{i,2});
            S{1}=acc(1:2,:)'/60;
            S{2}=dec(1:2,:)'/60;

            plotters(i)=fhrplot({FHR,FHRraw,baseline},TOCO,[analyses{i,1} ' - ' filename],S);

        end
        addlistener(plotters(i),'ChangeTime',@(src,evtdat) changeTime(i,plotters));
            
    end
end
function changeTime(n,Ps)
    for i=[1:n-1 n+1:length(Ps)]
        Ps(i).Time=Ps(n).Time;
        Ps(i).redraw()
    end
end

function Select(filename,BL,Accidents)
    g=load('expertAnalyses');
    ff=find(['\' filename]=='\'| ['\' filename]=='/',1,'last');
    s=find(strcmpi(filename(ff:end), {g.data(:).filename} ));
    if(isempty(s))
        s=size(g.Selection,1)+1;
    else
        g.data(s).trainingData=-1;
    end
    %s=35;
    a=dir(filename);
    d=(a.bytes-4)/6;
    baseline=linearinterpolation(BL(1,:)*240,BL(2,:),1:d);
   
    g.data(s)=struct('filename',filename(ff:end),'expertPts',BL,'baseline',baseline,'trainingData',g.data(s).trainingData,'accelerations',Accidents{1}, 'decelerations',Accidents{2}, 'overshoots', Accidents{3},'unreliableSignal', Accidents{4},'notToAnalyse',Accidents{5});
    save('expertAnalyses','-struct','g');
end