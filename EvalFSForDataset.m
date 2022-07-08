% Statistic evaluation of the FS model in the specific annotated dataset
%
% USAGE
%     EvalFSForDataset(model,stage,dataset,censorMHR)
%
% INPUT
%     model  : sensor to evaluate. Possible values: FSDop, FSScalp, FSMHR, or a
%              mat file with a cell table (n x 2) named PDop, PScalp OR PMHR. first
%              column the filename, second column the probabilities for
%              each sample of being an FS for the respective heart rate (F
%              or S)
%     stage  : Indicating if the evaluation should be done on 1: first
%              stage of delivery, 2: 2nd stage of delivery , 
%              -1 : Antepartum recordings, or 0: All recordings
%     dataset:  DopMHRTrain  (for Dopppler and MHR train dataset)
%               DopMHRVal    (for Dopppler and MHR Validation dataset)
%               ScalpTrain   (for Scalp ECG train dataset)
%               ScalpVal     (for Scalp ECG validation dataset)
%               DopMHRTestCP and DopMHRTestDbS ScalpTest Won't work as the labels are not provided
%     censorMHR true/false, if true the MHR is hidden to the model which
%              have to predict FS without it.
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

function [stats,results]=EvalFSForDataset(model,stage,dataset,censorMHR)
load('FSdataset/expertAnnotations','DB')
pernum=1;

stats.dataset=dataset;
stats.stage=stage;
stats.model=model;
stats.MHRcensored=censorMHR;

stats.MatchingPeriods=0;
stats.SigLen=0;
stats.NotLoss=0;
stats.Annotated=0;
stats.Accuracy=0;
stats.Sensibility=0;
stats.Specificity=0;
stats.PPV=0;
stats.PPN=0;
stats.AUC=0;
stats.CrossEntropy=0;
stats.WeightedAccuracy=0;
stats.WeightedCrossEntropy=0;
sumWeights=0;
AUCdata=zeros(2,0);
Contingence=zeros(2,2); % TN FN;FP TP
WeightedContingence=zeros(2,2);

if ~strcmp(model,'FSDop') && ~strcmp(model,'FSScalp') && ~strcmp(model,'FSMHR')
    DBm=load(model);
end
for i=1:length(DB)
    if strcmp(DB(i).dataset,dataset)
        if stage==0 || ...
                ( stage==-1 &&  strcmp(DB(i).Stage2_Start, 'Antepartum')) || ...
                ( stage==1 && (isempty(DB(i).Stage2_Start) || (isnumeric(DB(i).Stage2_Start) && DB(i).Stage2_Start>1) )) || ...
                ( stage==2 && ~isempty(DB(i).Stage2_Start) && isnumeric(DB(i).Stage2_Start) )

            d=dir(['FSDataset/**/' DB(i).filename]);
            [Sigs{2},Sigs{3},Sigs{1}]=fhropen([d.folder '/' d.name]);
            if censorMHR,Sigs{1}(:)=0;end
            isStage2=zeros(size(Sigs{2}));
            if ~isempty(DB(i).Stage2_Start) && isnumeric(DB(i).Stage2_Start)
                isStage2(DB(i).Stage2_Start:end)=1;
            end

            if strcmp(model,'FSMHR')
                [~,Prob]=FalseSigDetectDopMHR(Sigs{2},Sigs{1},isStage2);
                sig=1;
            elseif strcmp(model,'FSDop')
                Prob=FalseSigDetectDopMHR(Sigs{2},Sigs{1},isStage2);
                sig=2;
            elseif strcmp(model,'FSScalp')
                Prob=FalseSigDetectScalp(Sigs{3},isStage2);
                sig=3;
            else
                if isfield(DBm,'PDop')
                    matching=find(strcmp(DBm.PDop(:,1),DB(i).filename));
                    Prob=DBm.PDop{matching,2}; %#ok<*FNDSB> 
                    sig=2;
                elseif isfield(DBm,'PScalp')
                    matching=find(strcmp(DBm.PScalp(:,1),DB(i).filename));
                    Prob=DBm.PScalp{matching,2};
                    sig=3;
                elseif isfield(DBm,'PMHR')
                    matching=find(strcmp(DBm.PMHR(:,1),DB(i).filename));
                    Prob=DBm.PMHR{matching,2};
                    sig=1;
                end
            end


            TrueSel={DB(i).TS_MHR,DB(i).TS_FHR};TrueSel=TrueSel{min(sig,2)};
            FalseSel={DB(i).FS_MHR,DB(i).FS_FHR};FalseSel=FalseSel{min(sig,2)};
            TF=[TrueSel;FalseSel];
            if ~isempty(TF) && any(TF(:,2)>0 & TF(:,1)<length(Sigs{1}))
                s=1;
                e=length(Sigs{1});

                if stage==2  && ~isempty(DB(i).Stage2_Start) && isnumeric(DB(i).Stage2_Start)
                    Prob=Prob(DB(i).Stage2_Start:end);
                    TrueSel=TrueSel-DB(i).Stage2_Start+1;
                    FalseSel=FalseSel-DB(i).Stage2_Start+1;
                    s=round(DB(i).Stage2_Start);

                elseif stage==1 && ~isempty(DB(i).Stage2_Start) && isnumeric(DB(i).Stage2_Start)
                    Prob=Prob(1:round(DB(i).Stage2_Start)-1);
                    e=round(DB(i).Stage2_Start)-1;
                end
                S=Sigs{sig}(s:e);

                statper=evalPeriod(Prob,S,TrueSel,FalseSel);
                statper.dataNum=i;
                results(pernum)=statper; %#ok<AGROW>
                pernum=pernum+1;

                stats.SigLen=stats.SigLen+statper.SigLen;
                stats.NotLoss=stats.NotLoss+statper.NotLoss;
                stats.Annotated=stats.Annotated+statper.Annotated;
                stats.CrossEntropy=stats.CrossEntropy+statper.CrossEntropy*statper.Annotated;
                stats.WeightedCrossEntropy=stats.WeightedCrossEntropy+statper.WeightedCrossEntropy*statper.sumWeights;
                sumWeights=sumWeights+statper.sumWeights;
                Contingence=Contingence+statper.Contingence*statper.Annotated;
                WeightedContingence=WeightedContingence+statper.WeightedContingence*statper.sumWeights;
                AUCdata=[AUCdata statper.AUCdata]; %#ok<AGROW> 
            end
        end
    end
end
[~,~,~,stats.AUC]=perfcurve(AUCdata(1,:),AUCdata(2,:),1);
Contingence=Contingence/stats.Annotated;
stats.Sensibility=Contingence(2,2)/sum(Contingence(:,2));
stats.Specificity=Contingence(1,1)/sum(Contingence(:,1));
stats.PPV=Contingence(2,2)/sum(Contingence(2,:));
stats.PPN=Contingence(1,1)/sum(Contingence(1,:));
stats.Accuracy=sum(diag(Contingence));

WeightedContingence=WeightedContingence/sumWeights;
stats.WeightedAccuracy=sum(diag(WeightedContingence));

stats.CrossEntropy=stats.CrossEntropy/stats.Annotated;
stats.WeightedCrossEntropy=stats.WeightedCrossEntropy/sumWeights;


stats.MatchingPeriods=pernum-1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stats=evalPeriod(P,S,T,F)
P(S==0)=1;
stats.SigLen=length(S);
stats.NotLoss=sum(S>0);
Truth=.5*ones(size(S));
for i=1:size(T,1)
    Truth(max(1,round(T(i,1))):min(end,round(T(i,2))))=1;
end
for i=1:size(F,1)
    Truth(max(1,round(F(i,1))):min(end,round(F(i,2))))=0;
end


Truth(S==0)=.5;

Weights=ones(size(S))*sqrt(length(S)/sum(Truth~=.5));
Weights(Truth==.5)=0;
Annoted=(Truth~=.5);
stats.Annotated=sum(Annoted);
stats.sumWeights=sum(Weights);
if stats.Annotated>0
    stats.Contingence=[sum(Truth==1&P<0.5) sum(Truth==0&P<0.5);sum(Truth==1&P>=0.5) sum(Truth==0&P>=0.5)]/stats.Annotated;

    stats.WeightedContingence=[sum((Truth==1&P<0.5).*Weights) sum((Truth==0&P<0.5).*Weights);...
        sum((Truth==1&P>=0.5).*Weights) sum((Truth==0&P>=0.5).*Weights)]/stats.sumWeights;
    P=max(1.e-5,min(1-1.e-6,P));
    stats.CrossEntropy=-sum((log(P).*(Truth==0)+log(1-P).*(Truth==1)).*Annoted)/stats.Annotated;
    stats.WeightedCrossEntropy=-sum((log(P).*(Truth==0)+log(1-P).*(Truth==1)).*Weights)/stats.sumWeights;

    stats.AUCdata=[Truth(Weights>0);1-P(Weights>0)];
else
    stats.Contingence=zeros(2,2);
    stats.WeightedContingence=zeros(2,2);
    stats.CrossEntropy=0;
    stats.WeightedCrossEntropy=0;
    stats.AUCdata=zeros(2,0);
end
if any(Truth(Weights>0)==1) && any(Truth(Weights>0)==0)
    [~,~,~,stats.AUC]=perfcurve(Truth(Weights>0),1-P(Weights>0),1);
else
    stats.AUC=0;
end
stats.accuracy=sum(diag(stats.Contingence));
stats.WeightedAccuracy=sum(diag(stats.WeightedContingence));

stats.P=P;
stats.S=S;
stats.Y=Truth;

end