% Compute indices representing discordances between two analysis
%
% 
% USAGE
%    stats=statscompare(FHR,LDB1,LDB2,acc1,acc2,overshoots)
%         
% INPUT
%    acc*{1} acceleration acc*{2} deceleration
% OUTPUT
%    stats a struct
%       with all indexes :
%    MADI                   Morphological Analysis discordance index 
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


function stats=statscompare(FHR,LDB1,LDB2,acc1,acc2,overshoots)
    stats=struct();
    
    FHRu=FHR(FHR~=0);
    L1=LDB1(FHR~=0);
    L2=LDB2(FHR~=0);
    acc1{1}=filteraccident(acc1{1},FHR);
    acc2{1}=filteraccident(acc2{1},FHR);
    acc1{2}=filteraccident(acc1{2},FHR);
    acc2{2}=filteraccident(acc2{2},FHR);    
    overshoots=filteraccident(overshoots,FHR); 
    
    D1=conv((L1-FHRu).^2,ones(1,240)/240);D1=sqrt(D1(241:end-239))+3;
    D2=conv((L2-FHRu).^2,ones(1,240)/240);D2=sqrt(D2(241:end-239))+3;
    D=(L1(121:end-120)-L2(121:end-120)).^2;
    stats.MADI=mean(D./(D1.*D2+D));
    
    stats.RMSD_bpm=sqrt(mean((L1-L2).^2)); %Root Mean Square Difference
    
    stats.Diff_Over_15_bpm_prct=mean(abs(L1-L2)>15)*100;
    
    sL1=std(L1);sL2=std(L2);
    stats.Index_Agreement=1-sum((L1-L2).^2)/sum((abs(L1-mean(L1))+abs(L2-mean(L2))).^2);
    
    %Decelerations :
    [T,trueMatch,Qd]=accMatch(acc1{2},acc2{2},zeros(0,2));     
    stats.Dec_Match=T(1,1);
    stats.Dec_Only_1=T(1,3)+T(1,4);
    stats.Dec_Only_2=T(3,1)+T(4,1);
    stats.Dec_Doubled_On_1=T(1,4);
    stats.Dec_Doubled_On_2=T(4,1);
    stats.Dec_Only_1_Rate=(T(1,3)+T(1,4))/(T(1,1)+T(1,3)+T(1,4));
    if(T(1,1)+T(1,3)+T(1,4)==0)
        stats.Dec_Only_1_Rate=0;
    end
    stats.Dec_Only_2_Rate=(T(3,1)+T(4,1))/(T(1,1)+T(3,1)+T(4,1));
    if(T(1,1)+T(3,1)+T(4,1)==0)
        stats.Dec_Only_2_Rate=0;
    end    
    stats.Dec_Fmes=(1-stats.Dec_Only_1_Rate)*(1-stats.Dec_Only_2_Rate)/((1-stats.Dec_Only_1_Rate)+(1-stats.Dec_Only_2_Rate));
    
    stats.Dec_Start_RMSD_s=60*sqrt(mean((trueMatch(:,1)-trueMatch(:,3)).^2));
    stats.Dec_Start_Avg_2_M_1_s=-60*mean(trueMatch(:,1)-trueMatch(:,3));
    stats.Dec_Length_RMSD_s=60*sqrt(mean((trueMatch(:,2)-trueMatch(:,1)-trueMatch(:,4)+trueMatch(:,3)).^2));
    stats.Dec_Length_Avg_2_M_1_s=-60*mean( trueMatch(:,2)-trueMatch(:,1)-trueMatch(:,4)+trueMatch(:,3) );
    
    %Accelerations :
    [T,trueMatch,Qa]=accMatch(acc1{1},acc2{1},overshoots);
    stats.Acc_Match=T(1,1);
    stats.Acc_Only_1=T(1,3)+T(1,4);
    stats.Acc_Only_2=T(3,1)+T(4,1);
    stats.Acc_Doubled_On_1=T(1,4);
    stats.Acc_Doubled_On_2=T(4,1);
    stats.Acc_Only_1_Rate=(T(1,3)+T(1,4))/(T(1,1)+T(1,3)+T(1,4));
    if(T(1,1)+T(1,3)+T(1,4)==0)
        stats.Acc_Only_1_Rate=0;
    end
    stats.Acc_Only_2_Rate=(T(3,1)+T(4,1))/(T(1,1)+T(3,1)+T(4,1));
    if(T(1,1)+T(3,1)+T(4,1)==0)
        stats.Acc_Only_2_Rate=0;
    end
    stats.Acc_Fmes=(1-stats.Acc_Only_1_Rate)*(1-stats.Acc_Only_2_Rate)/((1-stats.Acc_Only_1_Rate)+(1-stats.Acc_Only_2_Rate));
    
    stats.Acc_Start_RMSD_s=60*sqrt(mean((trueMatch(:,1)-trueMatch(:,3)).^2));
    stats.Acc_Start_Avg_2_M_1_s=-60*mean(trueMatch(:,1)-trueMatch(:,3));
    stats.Acc_Length_RMSD_s=60*sqrt(mean((trueMatch(:,2)-trueMatch(:,1)-trueMatch(:,4)+trueMatch(:,3)).^2));
    stats.Acc_Length_Avg_2_M_1_s=-60*mean( trueMatch(:,2)-trueMatch(:,1)-trueMatch(:,4)+trueMatch(:,3) );
    
    %synthetic baseline inconsistency index
    FHRi=interpolFHR(FHR);
    cFHR1=cumsum(FHRi-LDB1)/240;cFHR2=cumsum(FHRi-LDB2)/240;
    AccArea1=[cFHR1(round(acc1{1}(:,2)*4*60))-cFHR1(round(acc1{1}(:,1)*4*60+1)) 0];
    AccArea2=[cFHR2(round(acc2{1}(:,2)*4*60))-cFHR2(round(acc2{1}(:,1)*4*60+1)) 0];
    DecArea1=[-cFHR1(round(acc1{2}(:,2)*4*60))+cFHR1(round(acc1{2}(:,1)*4*60+1)) 0];
    DecArea2=[-cFHR2(round(acc2{2}(:,2)*4*60))+cFHR2(round(acc2{2}(:,1)*4*60+1)) 0];
    
    d=0;dmax=0;dc=0;dmaxc=0;
    for i=1:size(Qa,1)
        for j=1:size(Qa,2)
            if(Qa(i,j)==1)
                d=d+(AccArea1(i)-AccArea2(j))^2;
                dmax=dmax+max([AccArea1(i) AccArea2(j)])^2;
                dc=dc+abs(AccArea1(i)-AccArea2(j));
                dmaxc=dmaxc+max([AccArea1(i) AccArea2(j)]);
            end
        end
    end
    stats.ASI_prct=sqrt(d)/sqrt(dmax)*100;
    %stats.ASIc_prct=dc/dmaxc*100;
    if(dmaxc==0)
        stats.ASI_prct=0;
        %stats.ASIc_prct=0;
    end    
    d=0;dmax=0;
    for i=1:size(Qd,1)
        for j=1:size(Qd,2)
            if(Qd(i,j)==1)
                d=d+(DecArea1(i)-DecArea2(j))^2;
                dmax=dmax+max([DecArea1(i) DecArea2(j)])^2;
                dc=dc+abs(DecArea1(i)-DecArea2(j));
                dmaxc=dmaxc+max([DecArea1(i) DecArea2(j)]);
            end
        end
    end
    stats.DSI_prct=sqrt(d)/sqrt(dmax)*100;
    %stats.DSIc_prct=dc/dmaxc*100;
    if(dmax==0)
        stats.DSI_prct=0;
        %stats.DSIc_prct=0;
    end
    stats.SI_prct=(stats.ASI_prct+2*stats.DSI_prct)/3;
    %stats.SIc_prct=(stats.ASIc_prct+2*stats.DSIc_prct)/3;
    
end

function [T,trueMatch,Q,counted1,counted2]=accMatch(a1,a2,e)
M1=cell(1,size(a1,1));

M2=cell(1,size(a2,1));
Q=zeros(size(a1,1),size(a2,1));

trueMatch=zeros(min([size(a1,1) size(a2,1)]),4);
n=1;

T=zeros(4);
counted1=zeros(1,size(a1,1));
counted2=zeros(1,size(a2,1));
for i=1:size(a1,1)
    M1{i}=find(a2(:,2)>a1(i,1)+5/60 & a2(:,1)+5/60<a1(i,2));
    Q(i,M1{i})=1;
    if(isempty(M1{i}))
        if any(abs(a1(i,1)-e(:,1))<15/60 & abs(a1(i,2)-e(:,2))<15/60)
            T(1,2)=T(1,2)+1;
            counted1(i)=-2;
        else
            T(1,3)=T(1,3)+1;
            counted1(i)=-1;
        end
    end
end
for i=1:size(a2,1)
    M2{i}=find(a1(:,2)>a2(i,1)+5/60 & a1(:,1)+5/60<a2(i,2));
    Q(M2{i},i)=1;
    if(isempty(M2{i}))
        if any(abs(a2(i,1)-e(:,1))<15/60 & abs(a2(i,2)-e(:,2))<15/60)
            T(2,1)=T(2,1)+1;
            counted2(i)=-2;
        else
            T(3,1)=T(3,1)+1;
            counted2(i)=-1;
        end
    elseif(length(M2{i})==1 && length(M1{M2{i}})==1 )
        T(1,1)=T(1,1)+1;
        counted2(i)=M2{i};
        counted1(M2{i})=i;
        trueMatch(n,:)=[a1(M2{i},1:2) a2(i,1:2)];
        n=n+1;
    end
end
trueMatch=trueMatch(1:n-1,:);

while ( any(counted1==0) || any(counted2==0) )
    for i=1:size(a1,1)
        if( isempty(M1{i}) && counted1(i)==0 )
            T(1,4)=T(1,4)+1;
            counted1(i)=-3;
        elseif( length(M1{i})==1  && counted1(i)==0 )
            for k=M2{M1{i}}'
                if k~=i
                    M1{k}=M1{k}(M1{k}~=M1{i});
                end
            end
            M2{M1{i}}=i; % If M2 as several matching
            if(counted2(M1{i})~=0)
                disp('Erreur 1');
            end
            T(1,1)=T(1,1)+1;
            counted1(i)=M1{i};
            counted2(M1{i})=i;
        end
    end
    for i=1:size(a2,1)
        if( isempty(M2{i}) && counted2(i)==0 )
            T(4,1)=T(4,1)+1;
            counted2(i)=-3;
        elseif( length(M2{i})==1  && counted2(i)==0 )
            for k=M1{M2{i}}'
                if k~=i
                    M2{k}=M2{k}(M2{k}~=M2{i});
                end
            end
            M1{M2{i}}=i; % If M2 as several matching
            
            if(counted1(M2{i})~=0)
                disp('Erreur 1');
            end
            T(1,1)=T(1,1)+1;
            counted2(i)=M2{i};
            counted1(M2{i})=i;
        end
    end
    
end

Q=Q(counted1~=-2, counted2~=-2);
Q(end+1,:)=(sum(Q)==0);
Q(:,end+1)=(sum(Q,2)==0);


end


function acc=filteraccident(acc,FHR)
keep=ones(1,size(acc,1));    
for i=1:size(acc,1)
    s=FHR( round(acc(i,1)*240+1):min([end round(acc(i,2)*240)]) )==0;
    if(sum(s)/length(s)>.33333)
        keep(i)=0;
    end
end
acc=acc(keep==1,:);
end