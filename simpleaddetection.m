% TRivial method for Acceleration/Deceleration detection
%
% USAGE
%    [acc,dec,falseacc,falsedec]=simpleaddetection(fhr,baseline)
%         
%
% INPUT
%     FHR       : Fetal Heart Rate sampled at 4Hz
%     baseline  : baseline sampled at 4Hz
%
% OUTPUT
%     acc       : 3xn for each n acceleration [start end max]
%     dec       : 3xn for each n deceleraiton [start end max]
%     falseacc  : Same for acceleration with only 5 bpm high
%     falsedec  : Same for deceleration with only 5 bpm high
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



function [acc,dec,falseacc,falsedec]=simpleaddetection(fhr,baseline)

acc=detectaccident(fhr-baseline,15);
dec=detectaccident(baseline-fhr,15);
falseacc=minusint(acc,detectaccident(fhr-baseline,5));
falsedec=minusint(dec,detectaccident(baseline-fhr,5));
end

function accidentsample=detectaccident(sig,thre)

peaks=find(sig>thre);
accidentsample=zeros(3,0);
while ~isempty(peaks)
    dacc=find(sig(1:peaks)<0,1,'last');
    if isempty(dacc)
        dacc=1;
    end
    facc=find(sig(dacc+1:end)<0,1,'first')+dacc;
    if isempty(facc)
        facc=length(sig);
    end
    [~,macc]=max(sig(dacc:facc));
    macc=macc+dacc-1;
    if facc-dacc>15*4
        accidentsample=[accidentsample [dacc;facc;macc]/4];
    end   
    peaks=peaks(peaks>facc);
end
end

function f=minusint(a,f)

for i=1:size(a,2)
    n=find(f(1,:)>=a(1,i) &f(2,:)<=a(2,i));
    if ~isempty(n)
        f=f(:,[1:n-1 n+1:end]);
    end
end
    
end