% Method MG
% Re-coded by S. Boudet from
% MONGELLI, Max, DAWKINS, Robert, CHUNG, Tony, et al. 
% Computerised estimation of the baseline fetal heart rate in labour: the low frequency line. 
% BJOG: an international journal of obstetrics & gynaecology, 1997, vol. 104, no 10, p. 1128-1133.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aammongelli(FHR)
%         Mongelli's method standard simple method for acceleration/deceleration detection 
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
function [baseline,acc,dec,falseacc,falsedec]=aammongelli(FHR)
sFHR=avgsubsamp(FHR,8);

Arrays=zeros(1,length(sFHR));
baseline=zeros(1,length(sFHR));
level=mod(sFHR(1:6*30),8);
baseline(1:3*30)=level;
Arrays(1:6*30)=-1;
Arrays(abs(sFHR(1:6*30)-level)<=8)=1;
% TO make the link with the paper "Main array" corresponds to points
% with Arrays==1 and "secondary array" corresponds to points with Arrays==-1

for winstart=1:length(sFHR)-6*30
    level=baseline(winstart+3*30-1);
    if(abs(sFHR(winstart+6*30)-level)<=8)
        Arrays(winstart+6*30)=1;
        %Paper : "Whenever a value is added into the main array, the oldest
        % fetal  heart  rate  values in both  arrays  are  discarded."
        % This is strange: it is much simpler, more justified and more efficient to work
        % on a 6 min windows. We guess it is an error on the description of the method
        % since it generates strange results
    else
        Arrays(winstart+6*30)=-1;
    end
    Arrays(winstart)=0;
    
    if sum(Arrays)<0
        [y,ny]=mod(sFHR(winstart+(1:6*30)),8);
        
        if((ny>=3*30 || ny>=2*sum(Arrays==1) )&& abs(y-level)>8)
            %ny>=sum(abs(sFHR(Arrays~=0)-level)<=8)
            % Condition : || ny>=2*sum(Arrays==1) is added because the program bug if there is not point with Arrays==1.
            % Anyway the main array would have no sense if there is too few
            % samples on it.
            
            
            %This part is incoherent on the paper.
            % On a first side : "the alternate array is switched to become the dominant"
            % On the other side : "The dominant array therefore always contains fetal heart
            % rate values within  a narrow band (<=8 bpm)  of a modal
            % value." The alternate array contains value >level+8bpm and
            % and others <level +8 bpm so the new main array cannot be
            % exactly the old alternate array. => We guess a new main array
            % and alternate array are computed from the modal value.
            Arrays(:)=0;
            Arrays(winstart+(1:6*30))=-1+2*(abs(sFHR(winstart+(1:6*30))-y)<=8);
            %fprintf('Switch t: %f, level:%f, new level:%f\n',winstart/30+3,level,mean(sFHR(Arrays==1)));
        end
    end
    baseline(winstart+3*30)=mean(sFHR(Arrays==1));
    
end
baseline(winstart+3*30:end)=baseline(winstart+3*30);
baseline=avgwin(baseline,6*30);
baseline=resamp(baseline,8,length(FHR));
[acc,dec,falseacc,falsedec]=simpleaddetection(FHR,baseline);
end


function [y,ny]=mod(points,n)
density=zeros(1,250);
for i=1:length(points)
    density(round(points(i))+(-n:n))=density(round(points(i))+(-n:n))+1;
end
[ny,y]=max(density);
end


