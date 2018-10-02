% Replace zeros in FHR corresponding to signal missing by linear
% interpolattion.
%
% USAGE
%     [FHR,d,f]=interpolFHR(FHR)
% INPUT 
%     FHR    : the FHR
% OUTPUT
%     FHR    : the interpolated FHR
%     d      : First no-missing signal sample
%     f      : Last no-missing signal sample
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
function [FHR,d,f]=interpolFHR(FHR)
    n=find(FHR>0,1);
    FHR(1:n)=FHR(n);
    d=n;
    while ~isempty(n) && n<length(FHR)
        n=find(FHR(n:end)==0,1)+n-1;
        nf=find(FHR(n:end)>0,1)+n-1;
        if(~isempty(nf))
            FHR(n-1:nf)=linspace(FHR(n-1),FHR(nf),nf-n+2);
        end
        n=nf;
    end
    
    f=find(FHR>0,1,'last');
    FHR(f:end)=FHR(f);
end