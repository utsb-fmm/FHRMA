% Provides the probabilities of being Fasle signals for FHR Scalp ECG
% channel thanks to FSScalp model
%
% USAGE
%     [PFHR,FHRtrue]=FalseSigDetectDopMHR(FHR,isStage2)
% INPUT 
%     FHR    : the Fetal Heart Rate (4 Hz) (line matrix)
%     isStage2: A vector either 0 for 1st stage of delivery samples or 
%               1 for second stage sample 
% OUTPUT
%     PFHR   : Probability for FHR signal of being False Signal
%     FHRtrue: FHR where samples with PFHR>0.5 replaced by Missing Signal
%     FHR: FHR where holds are removed
% 

%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     FHR Morphological Analysis Toolbox  Copyright (C) 2022 Samuel Boudet, Faculté de Médecine et Maïeutique,
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
function [PFHR,FHRtrue,FHR]=FalseSigDetectScalp(FHR,isStage2)
FHR=removeholds(FHR);

I=[(FHR>0).*(FHR-120)/60;FHR>0;isStage2];
[PFHR]=FSScalp(I);
FHRtrue=zeros(size(FHR));FHRtrue(PFHR<0.5)=FHR(PFHR<.5);

end



function FHR2=removeholds(FHR2)
%Find 12 consecutive identic samples and set them to 0

n=find(FHR2(1:end-4)~=FHR2(2:end-3)&FHR2(2:end-3)>0&...
    FHR2(2:end-3)==FHR2(3:end-2)&...
    FHR2(2:end-3)==FHR2(4:end-1)&...
    FHR2(2:end-3)==FHR2(5:end));
for i=n
    h=find(FHR2(i+2:min(length(FHR2),i+30*240))~=FHR2(i+1),1,'first');
    if h>=12
        %
        %std of consecutive beats [1;2.5] most of time; assume 0.5
        %P Delta consecutive beat diffence <.25bpm  =.2
        %if FHR>120bpm 1 sample over 2 is different form previous value
        %12 same values (3s) = 6 coincidences of equal consecutives beats each time p=.2
        %In 1h E(false hold)==60*240*p^(N/2)~1 with p=.2 and N=12 (3s)
        FHR2(i+1:i+h)=0;
    end
end


n=find(FHR2(1:end-1)>0&FHR2(2:end)==0);
for i=n
    h=find(FHR2(i-1:-1:max(i-17,1))~=FHR2(i),1,'first');
    if h>4
        FHR2(i-h+2:i)=0;
    end
end

n=find(FHR2(2:end-1)>0&FHR2(3:end)==0&FHR2(1:end-2)==0)+1;
FHR2(n)=0;

n=find(FHR2(2:end-2)>0&FHR2(2:end-2)==FHR2(3:end-1)&FHR2(4:end)==0&FHR2(1:end-3)==0)+1;
for i=n;FHR2(i:i+1)=0;end

n=find(FHR2(2:end-3)>0&FHR2(2:end-3)==FHR2(3:end-2)&FHR2(2:end-3)==FHR2(4:end-1)&FHR2(5:end)==0&FHR2(1:end-4)==0)+1;
for i=n;FHR2(i:i+2)=0;end
end

function [PDop]=FSScalp(I) 
I=I';
M=load('FSScalp');
RevI=I(end:-1:1,:,:);
IforwBack=cat(3,I,RevI);

GRU1=GRU(IforwBack,M.GRU1{:});
L1=cat(3,cat(2,GRU1(:,:,1:end/2),GRU1(end:-1:1,:,end/2+1:end),I),...
    cat(2,GRU1(end:-1:1,:,1:end/2),GRU1(:,:,end/2+1:end),RevI) );
GRU2=GRU(L1,M.GRU2{:});
L2=cat(3,cat(2,GRU2(:,:,1:end/2),GRU2(end:-1:1,:,end/2+1:end),I),...
    cat(2,GRU2(end:-1:1,:,1:end/2),GRU2(:,:,end/2+1:end),RevI) );
GRU3=GRU(L2,M.GRU3{:});
CDop=cat(2,GRU3(:,:,1:end/2),GRU3(end:-1:1,:,end/2+1:end));
PDop=Dense(CDop,M.Dense1{:})';

end

function O=GRU(I,W,U,B)
% Re coding of a GRU layer
O=zeros(size(I,1),size(U,1),size(I,3));
n=size(W,2)/3;
Wx=pagemtimes(I,W);
B=repmat(B,1,1,size(I,3));

for i=1:size(I,1)
    if i>1
        Uh=pagemtimes(O(i-1,:,:),U);
    else
        Uh=zeros(1,size(W,2),size(I,3));
    end
    z=1./(   1+exp(-( Wx(i,1:n,:)    +Uh(1,1:n,:)    +B(1,1:n,:)    +B(2,1:n,:) ))   );
    r=1./(   1+exp(-( Wx(i,n+1:2*n,:)+Uh(1,n+1:2*n,:)+B(1,n+1:2*n,:)+B(2,n+1:2*n,:) ))   );
    h=tanh(   Wx(i,2*n+1:3*n,:)+( Uh(1,2*n+1:3*n,:)+B(2,2*n+1:3*n,:) ).*r +B(1,2*n+1:3*n,:)   );
    if i>1
        O(i,:,:)=z.*O(i-1,:,:)+(1-z).*h;
    else
        O(i,:,:)=(1-z).*h;
    end
end
end

function O=Dense(I,W,b)
%Re-coding of a Fully connected layer with sigmoid activation
O=1./(  1+exp(-pagemtimes(I,W)-b)  );
end