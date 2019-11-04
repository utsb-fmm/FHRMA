% Convert a files analysed with a specific method to Web format
%
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
load analyses/expertAnalyses.mat

for i=1:156
    if i<=66
        f='traindata/';
    else
        f='testdata/';
    end
    
    
    [FHR1,FHR2,TOCO,timestamp]=fhropen([f data(i).filename]);
    FHRr=FHR2; FHRr(FHRr==0)=FHR1(FHRr==0);
    [FHRi,FHR,TOCOi,d,f]=preprocess(FHR1,FHR2,TOCO,data(i).unreliableSignal);
    
    [ baseline,ac,dc ] = aamwmfb(FHRi);
    accelerations=ac(1:2,:)'/60;
    decelerations=dc(1:2,:)'/60;
    filename=['web/' data(i).filename(1:end-3) 'WMFB.fhra'];
    
    
    
    FHRm=zeros(1,length(FHR));
    fhrsave(filename,FHR,FHRm,TOCOi,round(timestamp+d/4),FHRm,FHRi,baseline)
    
    [ ~,ct,~ ] = aamwmfb(TOCOi*2);
    contractions=ct(1:2,:)'/60;
    US=data(i).unreliableSignal;
    NTA=data(i).notToAnalyse;
    
    T=cell(0,2);
    for j=1:size(accelerations,1)
        T{j,1}=round(accelerations(j,1)*240);
        T{j,2}=['$$ACC ' num2str(round(accelerations(j,2)*240-accelerations(j,1)*240))];
    end
    N=size(T,1);
    for j=1:size(decelerations,1)
        T{N+j,1}=round(decelerations(j,1)*240);
        T{N+j,2}=['$$DEC ' num2str(round(decelerations(j,2)*240-decelerations(j,1)*240))];
    end
    N=size(T,1);
    for j=1:size(contractions,1)
        T{N+j,1}=round(contractions(j,1)*240);
        T{N+j,2}=['$$CON ' num2str(round(contractions(j,2)*240-contractions(j,1)*240))];
    end
    N=size(T,1);
    for j=1:size(US,1)
        T{N+j,1}=round(US(j,1)*240);
        T{N+j,2}=['$$URS ' num2str(round(US(j,2)*240-US(j,1)*240))];
    end
    N=size(T,1);
    for j=1:size(NTA,1)
        T{N+j,1}=round(NTA(j,1)*240);
        T{N+j,2}=['$$NTA ' num2str(round(NTA(j,2)*240-NTA(j,1)*240))];
    end
    T=sortrows(T,1);
    fid=fopen([filename 'h'],'w');
    for j=1:size(T,1)
        fprintf(fid,'%07d %s\n',T{j,1},T{j,2});
    end
    fclose(fid);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i<=66)
        baseline=data(i).baseline;
        accelerations=data(i).accelerations;
        decelerations=data(i).decelerations;
        filename=['web/' data(i).filename(1:end-3) 'consensus.fhra'];
        fhrsave(filename,FHR,FHRm,TOCOi,round(timestamp+d/4),FHRm,FHRi,baseline)
        US=data(i).unreliableSignal;
        NTA=data(i).notToAnalyse;
        
        T=cell(0,2);
        for j=1:size(accelerations,1)
            T{j,1}=round(accelerations(j,1)*240);
            T{j,2}=['$$ACC ' num2str(round(accelerations(j,2)*240-accelerations(j,1)*240))];
        end
        N=size(T,1);
        for j=1:size(decelerations,1)
            T{N+j,1}=round(decelerations(j,1)*240);
            T{N+j,2}=['$$DEC ' num2str(round(decelerations(j,2)*240-decelerations(j,1)*240))];
        end
        N=size(T,1);
        for j=1:size(US,1)
            T{N+j,1}=round(US(j,1)*240);
            T{N+j,2}=['$$URS ' num2str(round(US(j,2)*240-US(j,1)*240))];
        end
        N=size(T,1);
        for j=1:size(NTA,1)
            T{N+j,1}=round(NTA(j,1)*240);
            T{N+j,2}=['$$NTA ' num2str(round(NTA(j,2)*240-NTA(j,1)*240))];
        end
        T=sortrows(T,1);
        fid=fopen([filename 'h'],'w');
        for j=1:size(T,1)
            fprintf(fid,'%07d %s\n',T{j,1},T{j,2});
        end
        fclose(fid);
        
        fprintf('.');
    end
end
