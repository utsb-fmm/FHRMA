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



function plotters=fhrmorphoopenandplot(filename,analyses,showExpertReference)
    g=load('FHRMAdataset/analyses/expertAnalyses');
    if nargin<=2,showExpertReference=false;end
    ff=find(['\' filename]=='\'| ['\' filename]=='/',1,'last');
    s=find(strcmpi(filename(ff:end), {g.data(:).filename} ));
    
    [FHR1,FHR2,MHR,TOCOorig]=fhropen(filename);
    [FHR,FHRraw,TOCO,d,f]=preprocess(FHR1,FHR2,TOCOorig);
    RefBLPoints=zeros(2,0);
    ExpertsPts=zeros(2,0);
    if ~isempty(s)
        [FHR,FHRraw2,TOCO]=preprocess(FHR1,FHR2,TOCOorig,g.data(s).unreliableSignal);
        FHR0=FHRraw2;FHR0(isnan(FHR0))=0; %Used for stats
        rjct=[g.data(s).notToAnalyse;[length(FHR)/240-1 length(FHR)/240]] ;
        for j=1:size(rjct,1)
            FHR0(round(rjct(j,1)*240+1):round(rjct(j,2)*240))=0;
        end

        ExpertAcc={g.data(s).accelerations, g.data(s).decelerations, g.data(s).overshoots, g.data(s).unreliableSignal, g.data(s).notToAnalyse};
        ExpertsPts=g.data(s).expertPts;
        if showExpertReference
            RefSelectionAcc=ExpertAcc;
            RefBLPoints=g.data(s).expertPts;
        else
            RefSelectionAcc={zeros(0,2), zeros(0,2), zeros(0,2), g.data(s).unreliableSignal, g.data(s).notToAnalyse};        
        end
    else
        FHR0=FHRraw;FHR0(isnan(FHR0))=0;
        ExpertAcc(1:5)={zeros(0,2)};
        RefSelectionAcc(1:5)={zeros(0,2)};
    end
    
    for i=1:size(analyses,1)
        if strcmp(analyses{i,2},'expert') 
            if ~isempty(s), RefBLPoints=g.data(s).expertPts;else, RefBLPoints=zeros(2,0);end
            plotters(i)=fhrplotMA({FHRraw,FHR,MHR},TOCO,[analyses{i,1} ' - ' filename],ExpertAcc,ExpertsPts); %#ok<*AGROW>
            addlistener(plotters(i),'Save',@(src,evtdat) Save(filename,plotters(i).RefBLPoints,plotters(i).Selections));
        else
            [baseline,acc,dec]=eval(analyses{i,2});
            S{1}=acc(1:2,:)'/60;
            S{2}=dec(1:2,:)'/60;            
            acc=round(acc*4);dec=round(dec*4);

            MASigs=[baseline;baseline;baseline];
            for j=1:size(acc,2)
                MASigs(3,acc(1,j)+1:acc(2,j))=FHR(acc(1,j)+1:acc(2,j));
            end
            for j=1:size(dec,2)
                MASigs(2,dec(1,j)+1:dec(2,j))=FHR(dec(1,j)+1:dec(2,j));
            end

            plotters(i)=fhrplotMA({FHRraw,FHR,MHR,MASigs},TOCO,[analyses{i,1} ' - ' filename],RefSelectionAcc,RefBLPoints);
            if ~isempty(s) && ~isempty(g.data(s).baseline)
                stats=statscompare(FHR0(1:f-d+1),g.data(s).baseline(1:f-d+1),baseline,RefSelectionAcc(1:2),S,g.data(s).overshoots);
                plotters(i).StatText=sprintf(' MADI:%0.2f%% \n SI:%0.2f%%\n\n RMSD baseline: %.2f bpm \n 15bpm difference rate:  %0.2f%% \n\n Deceleration :\n Sensitivity: %.3f\n PPV: %.3f\n F-measure: %.3f\n RMSD durations: %.2f s\n Mean diff duration: %.3f s \n\n Acceleration :\n Sensitivity: %.3f\n PPV: %.3f\n F-measure: %.3f\n RMSD durations: %.2f s\n Mean diff duration: %.3f s',...
                    stats.MADI*100,stats.SI_prct,stats.RMSD_bpm,stats.Diff_Over_15_bpm_prct,1-stats.Dec_Only_1_Rate,1-stats.Dec_Only_2_Rate,2/(1/(1-stats.Dec_Only_2_Rate)+1/(1-stats.Dec_Only_1_Rate)) ,stats.Dec_Length_RMSD_s,stats.Dec_Length_Avg_2_M_1_s,1-stats.Acc_Only_1_Rate,1-stats.Acc_Only_2_Rate,2/(1/(1-stats.Acc_Only_2_Rate)+1/(1-stats.Acc_Only_1_Rate)),stats.Acc_Length_RMSD_s,stats.Acc_Length_Avg_2_M_1_s);
            end
        end
    end

    for i=1:length(plotters)
        addlistener(plotters(i),'ChangeTime',@(src,evtdat) changeTime(i,plotters));
    end
end

function changeTime(n,Ps)
    for i=[1:n-1 n+1:length(Ps)]
        Ps(i).Time=Ps(n).Time;
        Ps(i).redraw()
    end
end

function Save(filename,BL,Accidents)
    g=load('FHRMAdataset/analyses/expertAnalyses.mat');
    ff=find(['\' filename]=='\'| ['\' filename]=='/',1,'last');
    s=find(strcmpi(filename(ff:end), {g.data(:).filename} ));
    if(isempty(s))
        s=size(g.Selection,1)+1;
    else
        g.data(s).trainingData=-1;
    end
    [FHR1]=fhropen(filename);
    baseline=linearinterpolation(BL(1,:)*240,BL(2,:),1:length(FHR1));
   
    g.data(s)=struct('filename',filename(ff:end),'expertPts',BL,'baseline',baseline,'trainingData',g.data(s).trainingData,'accelerations',Accidents{1}, 'decelerations',Accidents{2}, 'overshoots', Accidents{3},'unreliableSignal', Accidents{4},'notToAnalyse',Accidents{5});
    save('FHRMAdataset/analyses/expertAnalyses.mat','-struct','g');
end