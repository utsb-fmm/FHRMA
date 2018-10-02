% Method C
% Re-coded by S. Boudet from
% S. Cazares - Automated Identification of Abnormal Patterns in the
% Intrapartum Cardiotocogram - phD thesis 2002 Oxford University
%
% USAGE
%    [baseline,accelerations,decelerations]=aamcazares(FHR)
%         Cazares's method with its method for acceleration/deceleration detection
%    [baseline,accelerations,decelerations]=aamcazares(FHR,1)
%         Cazares's method with a standard simple method for acceleration/deceleration detection
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
%     FHR Morphological Analysis Toolbox Copyright (C) 2018 Samuel Boudet, Faculté de Médecine et Maïeutique,
%     samuel.boudet@gmail.com
%
%     This file is part of FHR Morphological Analysis Toolbox
%
%     FHR Morphological Analysis Toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     FHR Morphological Analysis Toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [baseline,acc,dec,falseacc,falsedec]=aamcazares(FHR,simpleacc)
    ssfhr=avgsubsamp(FHR,15);   
    % OCS Method
    baseline=avgwin(close(open(ssfhr,19),137),65);
    acc=zeros(3,0);
    falseacc=zeros(3,0);
    falsedec=zeros(3,0);
    
    dec=deceleration(ssfhr,baseline)*15/4;
    acc=deceleration(240-ssfhr,240-baseline)*15/4;
    
    baseline=resamp(baseline,15,length(FHR));
    if( exist('simpleacc','var') && simpleacc==1)
        [acc,dec]=simpleaddetection(FHR,baseline);
    end
    
    
end

function s2=erode(s1,w)
    s2=zeros(1,length(s1));
    for i=1:length(s1)
        win=max(1,i-w):min(length(s1),i+w);
        s2(i)=min(s1(win));
    end
end

function s2=dilate(s1,w)
    s2=zeros(1,length(s1));
    for i=1:length(s1)
        win=max(1,i-w):min(length(s1),i+w);
        s2(i)=max(s1(win));
    end
end

function s2=open(s1,w)
    s2=dilate(erode(s1,w),w);
end
function s2=close(s1,w)
    s2=erode(dilate(s1,w),w);
end



function deceleration=deceleration(FHR,baseline)
    vFHR=open(FHR-baseline,19);
    dvFHR=(vFHR(2:end)-vFHR(1:end-1))*16;%unit bpm/min
    
    %Setting thresholds at 5th and 95th percentiles of the
    %derivate other a sliding window
    Winsize=27*16;
    Step=16;
    n=length(dvFHR);
    thresholds5p=zeros(1,n);
    thresholds95p=zeros(1,n);
    winInputStart=[1:Step:n-Winsize n-Winsize+1];
    winInputEnd=winInputStart+Winsize-1;
    
    winOutputStart=[1 winInputStart(2:end)+(Winsize-Step)/2];
    winOutputEnd=[winInputStart(1:end-1)+(Winsize+Step)/2-1 n];
    for i=1:length(winInputStart)
        s=sort(dvFHR(winInputStart(i):winInputEnd(i)));
        thresholds5p(winOutputStart(i):winOutputEnd(i))=s(22);
        thresholds95p(winOutputStart(i):winOutputEnd(i))=s(410);
    end
    
    %Deceleration start = first value where derivate under the 5th percentiles
    %Deceleration end    = first value where  derivate over the 95th percentiles
    dipStartCandidat=find(dvFHR(2:end)<thresholds5p(2:end) & dvFHR(1:end-1)>=thresholds5p(2:end))+1;
    dipEndCandidat=find(dvFHR(2:end)>thresholds95p(2:end) & dvFHR(1:end-1)<=thresholds95p(2:end))+1;
    
    %Identifying couples Start/End which define decelerations :
    %A Start is a first start candidat after an End and a End is a last end
    %candidat before a start
    deceleration=zeros(3,0);
    i=1;
    while isempty(i) || i<=length(dipStartCandidat) 
        startDIP=dipStartCandidat(i);
        firstendnext=find(dipEndCandidat>startDIP,1,'first');
        if ~isempty(firstendnext)
            
            i=find(dipStartCandidat>dipEndCandidat(firstendnext),1,'first');
            if isempty(i)
                i=length(dipStartCandidat) +1;
                e=length(dipEndCandidat);
            else
                e=find(dipEndCandidat<dipStartCandidat(i),1,'last');
            end
            endDIP=dipEndCandidat(e);
            [~,minDIP]=min(FHR(startDIP:endDIP));
            minDIP=minDIP+startDIP-1;
            deceleration=[deceleration [startDIP;endDIP;minDIP]];
        else
            i=length(dipStartCandidat) +1;
        end
    end
    
    %Remove decelerations where area<6 beats
    C=cumsum(baseline-FHR)/16; %Primitive of the signal
    Areas=C(deceleration(2,:))-C(deceleration(1,:));
    deceleration=deceleration(:,Areas>=6);

end