% Standard Preprocessing of FHR signal
%   Merge FHR1 and FHR2
%   Remove empty signal at begining and end
%   Remove small part signal which are untrustable
%   Linear interpollation of missing part.
%
% USAGE
%    [FHRi,FHR,TOCOi,d,f]=preprocess(FHR1,FHR2,TOCO)
% OUTPUT
%       d     : first sample of non missing signal
%       f     : last sample of non missing signal
%       The following output signals correspond to sample d to f of FHR1-2
%       and TOCO
%       FHR   :the FHR signal with NAN on missing signal
%       FHRi  :the FHR signal with linearn interpollation on missing signal
%       TOCOi :TOCO
%
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
function [FHRi,FHR,TOCO,d,f]=preprocess(FHR1,FHR2,TOCO,unreliableSignal)

FHR=max([FHR1;FHR2]);
if nargin>=4
    for j=1:size(unreliableSignal,1)
        FHR(round(unreliableSignal(j,1)*240+1):round(unreliableSignal(j,2)*240))=0;
    end
end
FHR=removesmallpart(FHR);
d=find(FHR>0,1);


[FHRi,d,f]=interpolFHR(FHR);
%FHRi=FHRi(d:f);TOCOi=TOCO(d:f);FHR=FHR(d:f);
FHR(FHR==0)=NaN;