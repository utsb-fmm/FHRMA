% Method A
% Re-coded by S. Boudet from
% D. Ayres-de Campos, J. Bernardes, A. Garrido, J. Marques-de-Sá, L. Pereira-Leite, 
% SisPorto 2.0: a program for automated analysis of cardiotocograms, 
% J. Matern. Fetal Med. 9 (2000) 311–318.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aamayres(FHR)
%         Ayres' method with standard simple method for acceleration/deceleration detection
%         
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
function [baseline,acc,dec,falseacc,falsedec]=aamayres(fhr)

winsize=10;
winstep=5;
srate=240;
baseline=zeros(size(fhr));
BL=baselinegraph(fhr(1:winsize*srate));
baseline(1:(winsize+winstep)*srate/2)=BL;
for i=1:ceil((length(fhr)/srate-winsize)/winstep)
    BL=baselinegraph(fhr(1+(i-1)*winstep*srate:(winsize+(i-1)*winstep)*srate));
    baseline((i-1)*winstep*srate+(winsize-winstep)*srate/2+(1:winstep*srate))=BL;
end
BL=baselinegraph(fhr(end-winsize*srate+1:end));
baseline(end-(winsize+winstep)*srate/2+1:end)=BL;
[acc,dec,falseacc,falsedec]=simpleaddetection(fhr,baseline);
end

function BL=baselinegraph(fhr)
[aSTV]=STV(fhr);
fhr1=avgwin(fhr,5);
[h,f]=histog(fhr1);

BL=f(1);
if BL>=110
    if BL>152
        for i=2:length(f)
            if(f(i)<BL && f(i)>110 && h(i)>1.6*aSTV*h(1)) 
                BL=f(i);
            end
        end
    else
        if aSTV<0.2
            F=4;
        elseif aSTV<0.3
            F=2;
        elseif aSTV<0.4
            F=1;
        elseif aSTV<0.6
            F=.5;
        else
            F=1;
        end
        
        for i=2:length(f)
            if(f(i)<BL && f(i)>110 && h(i)>F*aSTV*h(1)) % Clair problème
                BL=f(i);
            end
        end
    end
else
    for i=2:length(f)
        if(f(i)>110 && h(i)>(1-aSTV)/3*h(1)) 
            BL=f(i);
            return;
        end
    end
    for i=1:length(f)
        if(f(i)<BL && h(i)>aSTV*h(1)) 
            BL=f(i);
        end
    end
end

end

function [aSTV]=STV(fhr)
    aSTV=sum(abs(fhr(2:end)-fhr(1:end-1))<=1)/(length(fhr)-1);
end

function [h,f]=histog(fhr)
fhr=round(fhr);
histo=zeros(1,240);
for i=1:length(fhr)
    histo(fhr(i))=histo(fhr(i))+1;
end
histo=histo/sum(histo);
[h,f]=sort(histo,'descend');
f=f(1:50);h=h(1:50);
f=f(h>=0.05);
h=h(h>=0.05);

if(isempty(f))
    [h,f]=max(histo);
end

end