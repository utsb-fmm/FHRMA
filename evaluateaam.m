% Compute indices representing discordances between two analysis
%
% 
% USAGE
%    fileresults=evaluateaam(command,expertbase,folder)
% EXAMPLE
%    load analyses/expertAnalyses.mat   
%    fileresults=evaluateaam('aamlu(FHR)',data([data(:).trainingData]==1),'traindata/')
%       ->Evaluate the Method Lu on the training dataset
%       The main criteria is then MADI=median([fileresults(:).MADI]);
%    fileresults=evaluateaam('analyses/L_std.mat',data([data(:).trainingData]==1),'traindata/')
%       ->Same since Lu results are saves in L_std.mat         
% INPUT
%    command    EITHER the command to get the analysis
%                           ex 'aamlu(FHR)' ; 
%                             FHR is the interpollated FHR;
%                             FHR0 is the raw FHR with 0 on missing signal
%                             TOCO is the TOC
%                             Each command should return
%                             [baseline,acceleration,deceleration]
%               OR the mat file with a saved analyses (ex. L_std.mat)
%
%   expertbase  database of expert analysis, syntax used on expertAnalyses.mat   
%   folder      fhr file forlder

% OUTPUT
%   fileresults One line for each FHR file with the following indices
%    MADI                   Morphological Analysis Discordance index 
%    RMSD_bpm               Root Mean square difference of baseline
%    Diff_Over_15_bpm_prct  percent of time on which there is 15bpm baseline diffferences
%    Index_Agreement        Removed index
%    Dec_Match              Number of common deceleration
%    Dec_Only_1             Number of deceleration detected only by Analyse 1
%    Dec_Only_2             Number of deceleration detected only by Analyse 2
%    Dec_Doubled_On_1       Number of deceleration detected on Analyse 2 corresponding to two decelerations on Analyse 1
%    Dec_Doubled_On_2       Number of deceleration detected on Analyse 1 corresponding to two decelerations on Analyse 2 
%    Dec_Only_1_Rate        Rate of deceleration detected only by Analyse 1
%    Dec_Only_2_Rate        Rate of deceleration detected only by Analyse 2
%    Dec_Fmes               F-measure for deceleration
%    Dec_Start_RMSD_s       Root mean square difference for begining of deceleration (s)
%    Dec_Start_Avg_2_M_1_s  average difference for begining of deceleration (s)
%    Dec_Length_RMSD_s      Root mean square difference for length of deceleration (s)
%    Dec_Length_Avg_2_M_1_s average difference for length of deceleration (s)
%    Acc_Match              Number of common acceleration
%    Acc_Only_1             Number of acceleration detected only by Analyse 1
%    Acc_Only_2             Number of acceleration detected only by Analyse 2
%    Acc_Doubled_On_1       Number of acceleration detected on Analyse 2 corresponding to two accelerations on Analyse 1
%    Acc_Doubled_On_2       Number of acceleration detected on Analyse 1 corresponding to two accelerations on Analyse 2 
%    Acc_Only_1_Rate        Rate of acceleration detected only by Analyse 1
%    Acc_Only_2_Rate        Rate of acceleration detected only by Analyse 2
%    Acc_Fmes               F-measure for acceleration
%    Acc_Start_RMSD_s       Root mean square difference for begining of acceleration (s)
%    Acc_Start_Avg_2_M_1_s  average difference for begining of acceleration (s)
%    Acc_Length_RMSD_s      Root mean square difference for length of acceleration (s)
%    Acc_Length_Avg_2_M_1_s average difference for length of acceleration (s)
%    ASI_prct               Jezeski Synthetic Inconsistency for acceleration  
%    DSI_prct               Jezeski Synthetic Inconsistency for deceleration 
%    SI_prct                Jezeski Synthetic Inconsistency (global)
    

%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     FHR Morphological Analysis Toolbox  Copyright (C) 2018 Samuel Boudet, Faculté de Médecine et Maïeutique,
%     samuel.boudet@gmail.com
%
%     This file is part of "Fetal Heart Rate morphological analysis Toolbox"
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
function  fileresults=evaluateaam(command,expertbase,folder)
        if strcmp(command(end-3:end),'.mat')
            methodR=load(command,'data');
        end
    
    for i=1:length(expertbase)
        
        if nargin>=3
            filename=[folder expertbase(i).filename];
        else
            filename=expertbase(i).filename;
        end
        [FHR1,FHR2,TOCO]=fhropen(filename);
        [FHR,FHRraw,TOCO]=preprocess(FHR1,FHR2,TOCO,expertbase(i).unreliableSignal);
        FHR0=FHRraw;FHR0(isnan(FHR0))=0;

        
        
        
        if strcmp(command(end-3:end),'.mat')
            s=find(strcmp(expertbase(i).filename,{methodR.data(:).filename}));
            ADMethod={methodR.data(s).accelerations methodR.data(s).decelerations};
            baseline=methodR.data(s).baseline;
        else
            [baseline,accel,decel]=eval(command);
            ADMethod={accel(1:2,:)'/60;decel(1:2,:)'/60};
        end
        
        
        rjct=[expertbase(i).notToAnalyse;[length(FHR)/240-1 length(FHR)/240]] ;
        for j=1:size(rjct,1)
            FHR0(round(rjct(j,1)*240+1):round(rjct(j,2)*240))=0;
        end     
        ADExpert={expertbase(i).accelerations, expertbase(i).decelerations};
         

        %FHR0=FHR0(1:length(FHRraw)); %if rjct ends after, FHR0 size increased

        
        fileresults(i)=statscompare(FHR0,expertbase(i).baseline,baseline,ADExpert,ADMethod,expertbase(i).overshoots);
        fprintf('.');
    end
end