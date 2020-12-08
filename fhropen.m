% Open a binary *.fhr file or *.dat file (Physionet CTU-UHB database)
% 
% USAGE
%    [FHR1,FHR2,TOCO,timestamp,MHR]=fhropen(filename)
%
% INPUT
%     filename       : File location
%
% OUTPUT
%     FHR1           : First FHR signal (4Hz)
%     FHR2           : Second FHR signal (for twin or second sensor)
%     TOCO           : TOCO signal (4Hz)
%     timestamp      : Unix timestamp of the begining of recording
%     MHR            : Maternal Heart rate (if .fhrm)
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
function [FHR1,FHR2,TOCO,timestamp,MHR]=fhropen(filename)

f=fopen(filename,'r');

if(strcmp(filename(end-3:end),'.dat'))
    data=fread(f,[2,10000000],'uint16')/100;
    FHR1=data(1,:);
    FHR2=zeros(size(FHR1));
    TOCO=data(2,:);
    timestamp=0;
elseif(strcmp(filename(end-3:end),'.fhr') || strcmp(filename(end-3:end),'.rcf') )
    timestamp=fread(f,1,'uint32');

    data=fread(f,[3,10000000],'uint16');
    FHR1=data(1,:)/4;
    FHR2=data(2,:)/4;

    fseek(f,4,'bof');

    data=fread(f,[6,10000000],'uint8');
    TOCO=data(5,:)/2;
elseif(strcmp(filename(end-3:end),'.fhrm') || strcmp(filename(end-3:end),'.rcfm'))
    timestamp=fread(f,1,'uint32');

    data=fread(f,[4,10000000],'uint16');
    FHR1=data(1,:)/4;
    FHR2=data(2,:)/4;
    MHR=data(3,:)/4;

    fseek(f,4,'bof');

    data=fread(f,[6,10000000],'uint8');
    TOCO=data(5,:)/2;    
end

fclose(f);