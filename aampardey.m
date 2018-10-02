% Method MP
% Re-coded by S. Boudet from
% J. Pardey, M. Moulden, C.W.G. Redman, 
% A computer system for the numerical analysis of nonstress tests, 
% Am. J. Obstet. Gynecol. 186 (2002) 1095–1103.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aampardey(FHR)
%         Pardey's method with standard simple method for acceleration/deceleration detection
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
function [ baseline ,acc,dec,falseacc,falsedec] = aampardey( fhr )
hbi=round(60000./fhr);
subsamphbi=avgsubsamp(hbi,15);
baselinehbi=baselineloop(subsamphbi);
baselinefhr=60000./baselinehbi;
baseline=resamp(baselinefhr,15,length(fhr));

[acc,dec,falseacc,falsedec]=simpleaddetection(fhr,baseline);
end



function blhbi=baselineloop(hbi)
farfrombaseline=zeros(1,ceil(length(hbi)/32));

if(length(hbi)<1024)
    blhbi=baselinewindow(hbi,[-60 60]);
else
    blhbi=hbi;
    
    for i=1:ceil((length(hbi)-32)/32)
        deb=max([1-(i-1)*32 -495]);
        fin=min([528 length(hbi)-(i-1)*32]);
        bltmp=baselinewindow(hbi((i-1)*32+deb:(i-1)*32+fin),[-60 60]);
        d=min([32 length(blhbi)-(i-1)*32]);
        blhbi((i-1)*32+1:(i-1)*32+d)=bltmp(2-deb:-deb+d+1);
        if(all(hbi((i-1)*32+1:(i-1)*32+d)-bltmp(2-deb:-deb+d+1))<0)
            farfrombaseline(i)=1;
        elseif(all(hbi((i-1)*32+1:(i-1)*32+d)-bltmp(2-deb:-deb+d+1))>0)
            farfrombaseline(i)=-1;
        end
        if(i>=5&& all(farfrombaseline(i-4:i)>0))
            for j=i-4:i
                deb=max([1-(j-1)*32 -495]);
                fin=min([528 length(hbi)-(j-1)*32]);
                bltmp=baselinewindow(hbi((j-1)*32+deb:(j-1)*32+fin),[-60 120]);
                d=min([32 length(blhbi)-(j-1)*32]);
                blhbi((j-1)*32+1:(j-1)*32+d)=bltmp(2-deb:-deb+d+1);
            end
        end
        if(i>=5&& all(farfrombaseline(i-4:i)<0))
            for j=i-4:i
                deb=max([1-(j-1)*32 -495]);
                fin=min([528 length(hbi)-(j-1)*32]);
                bltmp=baselinewindow(hbi((j-1)*32+deb:(j-1)*32+fin),[-120 60]);
                d=min([32 length(blhbi)-(j-1)*32]);
                blhbi((j-1)*32+1:(j-1)*32+d)=bltmp(2-deb:-deb+d+1);
            end
        end        
    end
    
end

end
function bl=baselinewindow(hbi,margin)
p=baselineglobal(hbi);

hbi(hbi<p+margin(1))=p+margin(1);
hbi(hbi>p+margin(2))=p+margin(2);

bl=butterfilt(hbi,16,0,0.1,4,1);%Dawes : 0.1 ;Nyboe : 1.164
end
function p=baselineglobal(hbi)
%hbi=butterfilt(hbi,16,0,0.1,4,1);
hbi=round(hbi);
histogram=zeros(1,2000);
for i=1:length(hbi)
    if(~isnan(hbi(i)))
        histogram(hbi(i))=histogram(hbi(i))+1;
    end
end


histogram=conv([0.25 0.5 0.25],histogram);
histogram=histogram(2:end-1);
histogram=histogram/sum(histogram);
cshist=cumsum(histogram);
[~,mod]=max(histogram);

p=find(...
    cshist<=0.875 &...
    [0 0 0 0 0 (histogram(6:end-1)>histogram(7:end) &...
    histogram(6:end-1)>histogram(5:end-2) &...
    histogram(6:end-1)>histogram(4:end-3) &...
    histogram(6:end-1)>histogram(3:end-4) &...
    histogram(6:end-1)>histogram(2:end-5) &...
    histogram(6:end-1)>histogram(1:end-6) ...
    ) 0] & ...
    (abs((1:2000)-mod)<30 | histogram>0.05)... %Ligne ajouté de pardey
    ,1,'last');

if(isempty(p))
    p=mod;
end

end


