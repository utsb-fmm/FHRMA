% Method MD
% Re-coded by S. Boudet from
% MAEDA, Kazuo, UTSU, Masaji, NOGUCHI, Yasuaki, et al. 
% Central computerized automatic fetal heart rate diagnosis with a rapid and direct alarm system. 
% The Open Medical Devices Journal, 2012, vol. 4, no 1.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aammaeda(FHR)
%         Maeda's method with its method for acceleration/deceleration detection 
%         which is also the standard simple method for acceleration/deceleration detection
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
function [baseline,acc,dec,falseacc,falsedec]=aammaeda(FHR)

sFHR=avgsubsamp(FHR,8);
baseline=zeros(1,length(FHR));

for win=[0:150:length(sFHR)-151 length(sFHR)-150]
    
    bins=zeros(1,25);

    for i=1:150
        bins(ceil(sFHR(win+i)/10))=bins(ceil(sFHR(win+i)/10))+1;
    end
    [~,bestbins]=max(bins(1:20));
    
    baseline(win*8+1:win*8+1200)=mean(sFHR( sFHR<=bestbins*10 & sFHR>(bestbins-1)*10 ));

end


baseline(win*8+1201:length(FHR))=baseline(win*8+1200);


[acc,dec,falseacc,falsedec]=simpleaddetection(FHR,baseline);