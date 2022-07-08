% Provides the probabilities of being Fasle signals for FHR Doppler and MHR
% thanks to FSDop and FSMHR models
%
% USAGE
%     [PFHR,PMHR,FHRtrue,MHRtrue]=FalseSigDetectDopMHR(FHR,MHR,isStage2)
% INPUT 
%     FHR    : the Fetal Heart Rate (4 Hz) (line matrix)
%     MHR    : the Maternal Heart Rate
%     isStage2: A vector either 0 for 1st stage of delivery samples or 
%               1 for second stage sample 
% OUTPUT
%     PFHR   : Probability for FHR signal of being False Signal
%     PMHR   : Probability for MHR signal of being False Signal
%     FHRtrue: FHR where samples with PFHR>0.5 replaced by Missing Signal
%     MHRtrue: MHR where samples with MFHR>0.5 replaced by Missing Signal
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
%

function [PFHR,PMHR,FHRtrue,MHRtrue]=FalseSigDetectDopMHR(FHR,MHR,isStage2)

I=[(FHR>0).*(FHR-120)/60;FHR>0;(MHR>0).*(MHR-120)/60;MHR>0;isStage2];
[PFHR,PMHR]=FSDop(I);


FHRtrue=zeros(size(FHR));FHRtrue(PFHR<.5)=FHR(PFHR<.5);
MHRtrue=zeros(size(MHR));MHRtrue(PMHR<.5)=MHR(PMHR<.5);

end

function [PDop,PMat]=FSDop(I) 
% The Model architecture
I=I';
M=load('FSDop');
RevI=I(end:-1:1,:,:);
IforwBack=cat(3,I,RevI);
GRU1M=GRU(IforwBack,M.GRU1MHR{:});
IPMat=cat(2,GRU1M(:,:,1:end/2),GRU1M(end:-1:1,:,end/2+1:end));
PMat=Dense(IPMat,M.DensePmat{:});
RPMat=PMat(end:-1:1,:,:);
L0=cat(3,cat(2,PMat,I),cat(2,RPMat,RevI));
GRU1Dop=GRU(L0,M.GRU1{:});
L1=cat(3,cat(2,GRU1Dop(:,:,1:end/2),GRU1Dop(end:-1:1,:,end/2+1:end),PMat,I),...
    cat(2,GRU1Dop(end:-1:1,:,1:end/2),GRU1Dop(:,:,end/2+1:end),RPMat,RevI) );
GRU2Dop=GRU(L1,M.GRU2{:});
L2=cat(3,cat(2,GRU2Dop(:,:,1:end/2),GRU2Dop(end:-1:1,:,end/2+1:end),PMat,I),...
    cat(2,GRU2Dop(end:-1:1,:,1:end/2),GRU2Dop(:,:,end/2+1:end),RPMat,RevI) );
GRU3Dop=GRU(L2,M.GRU3{:});
CDop=cat(2,GRU3Dop(:,:,1:end/2),GRU3Dop(end:-1:1,:,end/2+1:end));
PDop=Dense(CDop,M.Dense1{:})';
PMat=PMat';
end


function O=GRU(I,W,U,B)
% Re coding of a GRU layer
O=zeros(size(I,1),size(U,1),size(I,3));
n=size(W,2)/3;
Wx=pagemtimes(I,W);%(Chaque dim1 et chaque dim3 de I sont indépendant)
%WxT=pagemtimes(WT,IT);%(Chaque dim2 et chaque dim3 de IT sont indépendant WxT(i,j,k)=WT(i,:,k)*IT(:,j,k) )
% Il est 10% plus rapide de rassembler les dimensions indépendantes
% (car A(n1,n2)*B(n2,n3) avec n3>>n1 et n2 est plus rapide que B'*A')
B=repmat(B,1,1,size(I,3));


for i=1:size(I,1)
    %Wx=I(i,:)*W;
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