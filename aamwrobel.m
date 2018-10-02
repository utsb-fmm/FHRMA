% Method W
% Re-coded by S. Boudet from
% J. Wróbel, K. Horoba, T. Pander, J. Je?ewski, R. Czaba?ski, 
% Improving fetal heart rate signal interpretation by application of myriad filtering, 
% Biocybern. Biomed. Eng. 33 (2013) 211–221.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aamwrobel(FHR)
%         Wrobel's method with standard simple method for acceleration/deceleration detection
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
function  [ baseline ,acc,dec,falseacc,falsedec] =aamwrobel(FHR)
    sFHR=avgsubsamp(FHR,10);
    baseline=zeros(1,length(sFHR));
    for i=1:length(sFHR)
        wins=max(1,i-160);
        wine=min(length(sFHR),i+160);
        baseline(i)=myriad(sFHR(wins:wine),0.25);
%         if(i==35*25)
%             myriad(sFHR(wins:wine),0.25);
%         end
    end
    baseline=resamp(baseline,10,length(FHR),1);
    
    [acc,dec,falseacc,falsedec]=simpleaddetection(FHR,baseline);
end



function m=myriad(x,k2)
    ma=round(max(x));
    mi=round(min(x));
    myrs=zeros(1,ma-mi+1);
    for i=0:ma-mi
        myrs(i+1)=sum(log(k2+(x-i-mi).^2));
    end
    [~,m]=min(myrs);
    m=m+mi;
end