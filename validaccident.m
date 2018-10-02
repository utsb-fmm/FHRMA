% Remove acceleration/deceleration which does not have duration and
% amplitude to be a true A/D

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


function [accidents,reject]=validaccident(accidents,varsig,durationthreshold,amplitudethreshold)

accidents=accidents(:,accidents(2,:)-accidents(1,:)>=durationthreshold);

t=zeros(1,size(accidents,2));

for i=1:size(accidents,2)
    t(i)=(  max( varsig(accidents(1,i)*4:accidents(2,i)*4) ) >=amplitudethreshold  );
end
reject=accidents(:,t==0);
accidents=accidents(:,t==1);