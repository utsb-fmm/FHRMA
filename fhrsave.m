% Save to a binary *.fhr file (Physionet CTU-UHB database)
% .MHR means there is maternal heart rate
% .fhra means there is a morphological analysis
%
% USAGE
%    fhrsave(filename,FHR1,FHR2,MHR,TOCO,timestamp,FHRi,baselineFHR)
%
% INPUT
%     filename       : File location
%
% OUTPUT
%     FHR1           : First FHR signal (4Hz)
%     FHR2           : Second FHR signal (for twin or second sensor)
%     TOCO           : TOCO signal (4Hz)
%     timestamp      : Unix timestamp of the begining of recording
%     MHR           : (optional) Maternal Heart Rate
%     FHRi           : (optional) Pr�processed FHR (interpollated)
%     baselineFHR    : (optional) Baseline FHR
%     baselineTOCO   : (optional) Baseline TOCO
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
function fhrsave(filename,FHR1,FHR2,MHR,TOCO,timestamp,infos,FHRi,baselineFHR)
f=fopen(filename,'w');
fwrite(f,timestamp,'uint32');

Q=zeros(1,length(FHR1));
if nargin>=7 && ~isempty(infos)
    Q=infos.Q1+2*infos.isECG1+4*infos.Q2+8*infos.isECG2+16*infos.Qm+32*infos.isTOCOMHR+64*infos.isIUP;
end

for i=1:length(FHR1)
    if ~isempty(MHR)
        fwrite(f,[FHR1(i);FHR2(i);MHR(i)]*4,'uint16');
    else
        fwrite(f,[FHR1(i);FHR2(i)]*4,'uint16');
    end
    fwrite(f,[TOCO(i)*2;Q(i)],'uint8');
    if nargin>=9
        fwrite(f,[FHRi(i);baselineFHR(i)]*4,'uint16');
    end
    
end
fclose(f);
end