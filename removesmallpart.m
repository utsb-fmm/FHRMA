% Remove aberrant sample of an FHR

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
function FHR=removesmallpart(FHR)
    FHR(FHR>220|FHR<50)=0;
    n=find(FHR(1:end-1)==0 & FHR(2:end)>0)+1;
    for i=1:length(n)
        f=find(FHR(n(i):end)==0,1,'first');
        if f<5*4
        	FHR(n(i):n(i)+f)=0;
        end
    end
    
    
    n=find(FHR(1:end-1)==0 & FHR(2:end)>0)+1;
    for i=1:length(n)
        f=find(FHR(n(i):end)==0,1,'first');
        if f<30*4
            
            lastvalid=find(FHR(1:n(i)-1)>0,1,'last');
            nextvalid=find(FHR(n(i)+f:end)>0,1,'first')+n(i)+f-1;
            
            try
                if(  (FHR(n(i))-FHR(lastvalid)<-25 && FHR(n(i)+f-2)-FHR(nextvalid)<-25 ))
                    FHR(n(i):n(i)+f)=0;
                end
                if(  (FHR(n(i))-FHR(lastvalid)>25 && FHR(n(i)+f-2)-FHR(nextvalid)>25 ))
                    FHR(n(i):n(i)+f)=0;
                end
            catch
            end
        end
    end

end