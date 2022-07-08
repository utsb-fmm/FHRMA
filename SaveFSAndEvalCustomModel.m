% This is a model function if you want to evaluate your own method.
% This example is for Dopppler FHR but you can change it for MHR or Scalp
% FHR
% It processes all the Train/validation/Test dataset and save the method
% analysis to MyFSDopMethodAnalysis.mat. Then it evaluate the analysis and
% Train and validation recording to check there is no encodage error.
%
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
function [stats,results]=SaveFSAndEvalCustomModel()

load('./FSdataset/expertAnnotations.mat','DB');

%Here we want to process files on dataset with name containing Dop
%Replace by Scalp if you want to process Scalp ECG recordings
DB=DB(contains({DB.dataset},'Dop'));

PDop=cell(length(DB),2); % Replace by PScalp or PMHR for those models

%Change to parfor if you have parrallel computing toolbox for faster
%computation
for i=1:length(DB) 
    % Open all files in the DB
    d=dir(['FSDataset/**/' DB(i).filename]);
    [FHRDop,FHRScalp,MHR,TOCO]=fhropen([d.folder '/' d.name]); %#ok<ASGLU> 
    
    %For FSDOP we need to prepare the binary signal to indicate which 
    % signal sample corresponds to second stage of delivery. This 
    % information is not stored in .fhrm files.

    isStage2=zeros(size(FHRDop));
    if ~isempty(DB(i).Stage2_Start) && isnumeric(DB(i).Stage2_Start)
        isStage2(DB(i).Stage2_Start:end)=1;
    end

    PDop{i,1}=d.name;
    PDop{i,2}=FalseSigDetectDopMHR(FHRDop,MHR,isStage2);
    fprintf('.')
    if rem(i,100)==0
        fprintf(' %d/%d\n',i,length(DB) );
    end
    
end
save('MyFSDopMethodAnalysis.mat','PDop')

%To be sure that the file proparation is correct we show evaluation on
%train/validation dataset; Result should be the same than EvalFSAllDatasets

AllStages=0;Stage1=1;Stage2=2;Antepartum=-1;
params={
    'MyFSDopMethodAnalysis.mat',AllStages,'DopMHRTrain';
    'MyFSDopMethodAnalysis.mat',AllStages,'DopMHRVal';%4
    'MyFSDopMethodAnalysis.mat',Antepartum,'DopMHRVal';
    'MyFSDopMethodAnalysis.mat',Stage1,'DopMHRVal';
    'MyFSDopMethodAnalysis.mat',Stage2,'DopMHRVal';
};


results=[];
stats=[];

for i=1:length(params)
    disp(params(i,:))
    [s,r]=EvalFSForDataset(params{i,1},params{i,2},params{i,3},false);
    stats=[stats s]; %#ok<AGROW> 
    results=[results r]; %#ok<AGROW> 
end

f=uifigure;uitable(f,'Data',struct2table(stats),'Unit','Normalized','Position',[0 0 1 1]);
end

