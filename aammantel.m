% Method MT
% Re-coded by S. Boudet from
% R. Mantel, H.P. van Geijn, F.J. Caron, J.M. Swartjes, E.E. van Woerden, H.W. Jongsma, 
% Computer analysis of antepartum fetal heart rate: 1. Baseline determination, 
% Int. J. Biomed. Comput. 25 (1990) 261–272. 
%
% USAGE
%    [baseline,accelerations,decelerations]=aammantel(FHR)
%         Mantel's method with its method for acceleration/deceleration detection
%    [baseline,accelerations,decelerations]=aammantel(FHR,1)
%         Mantel's method with a standard simple method for acceleration/deceleration detection
%
% INPUT
%     FHR       : Fetal Heart Rate sampled at 4Hz
%
% OUTPUT
%     baseline  : the baseline signal at 4Hz
%     accelerations : Table with begining and end of each accelerations in s in each column
%     decelerations : Table with begining and end of each decelerations in s in each column
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
function [baseline,accelerations,decelerations]=aammantel(FHR,FHR0,simpleacc)
   
    subsamphbi=avgsubsamp(round(60000./FHR),10); %Heart beat interval
    blsupport=baselineloop(subsamphbi);
    
    bl=filterpass(subsamphbi,blsupport);

    bl=trimRR(subsamphbi,bl,20,20);
    bl=filterpass(bl,blsupport);
    bl=trimRR(subsamphbi,bl,15,20);
    bl=filterpass(bl,blsupport);
    bl=trimRR(subsamphbi,bl,10,20);
    bl=filterpass(bl,blsupport);
    bl=trimRR(subsamphbi,bl,5,20);
    bl=filterpass(bl,blsupport);
    
    baselinefhr=60000./bl;
    
    fhr0subsamp=zeros(size(baselinefhr));
    for i=1:length(fhr0subsamp)
        FHR0w=FHR0((i-1)*10+1:i*10);
        fhr0subsamp(i)=mean(FHR0w(FHR0w>0));
    end
    fhr00subsamp=fhr0subsamp;fhr00subsamp(isnan(fhr00subsamp))=0;
    fhrsubsamp=interpolFHR(fhr00subsamp);
    accelerations=MantelAccident(fhrsubsamp-baselinefhr,isnan(fhr0subsamp));
    decelerations=MantelAccident(baselinefhr-fhrsubsamp,isnan(fhr0subsamp));
    
    baseline=resamp(baselinefhr,10,length(FHR));
    if( exist('simpleacc','var') && simpleacc==1)
        [accelerations,decelerations]=simpleaddetection(FHR,baseline);
    end
end


function acc=MantelAccident(s,signan)
    starts=find(s(1:end-1)<5 & s(2:end)>=5)+1;
    valids=zeros(size(starts));
    ends=valids;
    for i=1:length(starts)
        e=starts(i)+find(s(starts(i):end-1)>=5 & s(starts(i)+1:end)<5,1,'first')-1;
        if ~isempty(e)
            ends(i)=e;
        else
            ends(i)=length(s);
        end
        
        valids(i)=(ends(i)-starts(i)>=5) & any(s(starts(i):ends(i))>10);
    end
    
    valids=find(valids);
    %Try to enlarge accidents
    for i=1:length(valids)
        j=valids(i);
        if i==1;
            d=1;
        else
            d=valids(i-1)+1;
        end
        
        if i==length(valids)
            f=length(starts);
        else
            f=valids(i+1)-1;
        end
        D=j;F=j;
        for k=j-1:-1:d
            if( starts(k+1)-ends(k)<=5 && D==k+1 && all(s(ends(k):starts(k+1))>0) )
                D=k;
            end
        end
        for k=j+1:f
            if(starts(k)-ends(k-1)<=5 && F==k-1 && all(s(ends(k-1):starts(k))>0) )
                F=k;
            end
        end
        
        starts(j)=starts(D);
        ends(j)=ends(F);
        starts([D:j-1 j+1:F])=-10;
        ends([D:j-1 j+1:F])=-10;
    end
    
    %Phase II
    starts=starts(valids);
    ends=ends(valids);

    acc=zeros(0,2);
    for i=1:length(starts)
        starts(i)=starts(i)-1+find(~signan(starts(i):end),1,'first');
        if ~all(signan(ends(i):-1:1))
            ends(i)=ends(i)+1-find(~signan(ends(i):-1:1),1,'first');

            avgthre=max([1.00001 1/5*s(starts(i):ends(i))]);

            siglostStart=starts(i)+find(~signan(starts(i):ends(i)-1)& signan(starts(i)+1:ends(i)));
            siglostEnd=starts(i)-1+find(signan(starts(i):ends(i)-1)& ~signan(starts(i)+1:ends(i)));

            ncuts=find(siglostEnd-siglostStart+1>=avgthre);

            for j=1:length(ncuts)+1
                if j==1
                    D=starts(i);
                else
                    D=siglostEnd(ncuts(j-1))+1;
                end

                if j==length(ncuts)+1
                    F=ends(i);
                else
                    F=siglostStart(ncuts(j))-1;
                end

                acc(end+1,:)=[D F]; 
            end
        end
    end
    
    valids=zeros(1,size(acc,1));
    for i=1:length(valids)
        if(acc(i,2)-acc(i,1)>=5 && any(s(acc(i,1):acc(i,2))>10) )
            valids(i)=1;
        end
    end
    
    acc=acc(valids==1,:)'*2.5;
    
    
    
    
end



function blsupport=baselineloop(hbi)
blsupport=zeros(1,length(hbi));
%Suppose 64min window which is 64*24 =1536 samples
if(length(hbi)<1536)
    blsupport(:)=baselevelreference(hbi);
else
    for i=[0:1536:length(hbi)-1537 length(hbi)-1536]
        blsupport(i+1:i+1536)=baselevelreference(hbi(i+1:i+1536));
    end
end
end



function p=baselevelreference(hbi)
hbi=round(hbi);
histogram=zeros(1,300);
for i=1:length(hbi)
    if(~isnan(hbi(i)) && hbi(i)>300 && hbi(i)<=600)
        histogram(hbi(i)-300)=histogram(hbi(i)-300)+1;
    end
end

histogram=histogram/sum(histogram);
cshist=cumsum(histogram);

p=find(...
    cshist<=0.875 &...
    [0 0 0 0 0 (histogram(6:end-1)>histogram(7:end) &...
    histogram(6:end-1)>histogram(5:end-2) &...
    histogram(6:end-1)>histogram(4:end-3) &...
    histogram(6:end-1)>histogram(3:end-4) &...
    histogram(6:end-1)>histogram(2:end-5) &...
    histogram(6:end-1)>histogram(1:end-6) ...
    ) 0] ... 
    ,1,'last')+300;

if(isempty(p))
    [~,p]=max(histogram);
end

end


function B=filterpass(B,P)
    %Backward dummy
    B0=P(1);
    for i=length(B):-1:1
        if abs(B(i)-P(i))<=60
            B0=0.95*B0+0.05*B(i);
        end
    end
    
    %Forward pass
    if abs(B(1)-P(1))<=60
        B(i)=0.95*B0+0.05*B(i);
    else
        B(i)=B0;
    end
    for i=2:length(B)
        if abs(B(i)-P(i))<=60
            B(i)=0.95*B(i-1)+0.05*B(i);
        else
            B(i)=B(i-1);
        end
    end
    
    %Backward Pass
    for i=length(B)-1:-1:1
        B(i)=0.95*B(i+1)+0.05*B(i);
    end
end



function B=trimRR(A,B,U,L)
    blpoints=ones(1,length(A));
    A=60000./A;
    B=60000./B;
    i=1;
    while i<length(A)
        if A(i)>B(i)+U
            d=find(A(1:i-1)<B(1:i-1),1,'last')+1;
            if(isempty(d))
                d=1;
            end
            f=find(A(d:end)<B(d:end),1,'first')+d-2;
            blpoints(d:f)=0;
            i=f+1;
        elseif A(i)<B(i)-L
            d=find(A(1:i-1)>B(1:i-1),1,'last')+1;
            if(isempty(d))
                d=1;
            end
            f=find(A(d:end)>B(d:end),1,'first')+d-2;
            blpoints(d:f)=0;
            i=f+1;
        else
            i=i+1;
        end
    end
    B(blpoints==1)=A(blpoints==1);
    B=60000./B;
end
