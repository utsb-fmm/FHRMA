% Method Weigthed Median Filter Baseline (WMFB)
% Coded by S. Boudet from
% Boudet, S., Houzé de l’Aulnoit, A., Demailly, R., Peyrodie, L., Beuscart, R., Houzé de l’Aulnoit,D. -
% Fetal heart rate baseline computation with a weighted median filter.
% In : Computers in Biology and Medicine, In Press, 2019.
%
% USAGE
%    [baseline,accelerations,decelerations]=aamwmfb(FHR)
%
% INPUT
%     FHR       : Fetal Heart Rate sampled at 4Hz (with no signal loss)
%
% OUTPUT
%     baseline  : the baseline signal at 4Hz
%     accelerations : Table with begining and end of each accelerations in s in each column
%     decelerations : Table with begining and end of each decelerations in s in each column
%
%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     FHR Morphological Analysis Toolbox  Copyright (C) 2019 Samuel Boudet, Faculté de Médecine et Maïeutique,
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


function [baseline,accelerations,decelerations,falseAcc,falseDec]=aamwmfb(FHRi)

FHR1=butterfilt(FHRi,240,0,1,1);
FHR2=butterfilt(FHRi,240,0,2,1);
FHR4=butterfilt(FHRi,240,0,4,1);
FHR8=butterfilt(FHRi,240,0,8,1);
FHR16=butterfilt(FHRi,240,0,16,1);



fcut=[0 1 3 7];
fdat=cell(3,2);
for j=1:3
    fdat{j,1}=butterfilt(FHRi,240,fcut(j),fcut(j+1),1);
    fdat{j,2}=[0 fdat{j,1}(2:end)-fdat{j,1}(1:end-1)]*240;
end
t=[abs(fdat{1,2});...
    enveloppe(fdat{1,2},240,0,2*fcut(2));...
    enveloppe(fdat{2,2},240,0,2*fcut(3));...
    enveloppe(fdat{3,2},240,0,2*fcut(4))]';
Q=[-2.4744    0.0266    0.0413    0.0105    0.0036]';
P=1-exp(Q(1)+t*Q(2:end))./(1+exp(Q(1)+t*Q(2:end)));


distancecoef=[linspace(0,1,200) 1   linspace(1,0,200)];

[bl1,mp1]=medgliss(FHR2,distancecoef,P',24);
P2=P'.*exp(3.21-0.19*abs(FHR1-bl1))./(1+exp(3.21-0.19*abs(FHR1-bl1)) );
[bl2]=medgliss(FHR2,distancecoef.^2,P2,24,bl1,mp1,0.1);
P3=P'.*exp(2.5-0.19*abs(FHR4-bl2))./(1+exp(2.5-0.19*abs(FHR4-bl2)) );
[bl3]=medgliss(FHR4,distancecoef.^4,P3,24,bl2,mp1,0.1);
P4=P'.*exp(2-0.19*abs(FHR8-bl3))./(1+exp(2-0.19*abs(FHR8-bl3)) );
[bl4]=medgliss(FHR8,distancecoef.^8,P4,24,bl3,mp1,0.1);
P5=P'.*exp(1.5-0.19*abs(FHR16-bl4))./(1+exp(1.5-0.19*abs(FHR16-bl4)) );
[bl5]=medgliss(FHR16,distancecoef.^16,P5,24,bl4,mp1,0.1);
P6=P'.*exp(1-0.19*abs(FHR16-bl5))./(1+exp(1-0.19*abs(FHR16-bl5)) );
[bl6]=medgliss(FHR16,distancecoef.^16,P6,24,bl5,mp1,0.1);
baseline=bl6;


accelerations=accidentcandidat(FHR1,baseline,5);
accelerations=adjustduration(accelerations,FHRi-baseline)/4;
[accelerations,falseAcc]=validaccident(accelerations,FHRi-baseline,15,15);
decelerations=accidentcandidat(baseline,FHR1,5);
decelerations=adjustduration(decelerations,baseline-FHRi)/4;
[decelerations,falseDec]=validaccident(decelerations,baseline-FHRi,15,15);
end

%*************************************************************************
function startendlist=adjustduration(startendlist,s)
i=1;
while i<=size(startendlist,2)
    newend=startendlist(3,i)-1+find([s(startendlist(3,i):startendlist(2,i))<0 -1],1,'first');
    if (startendlist(2,i)-newend)>=15*4
        [m,imax]=max(s(newend+1:startendlist(2,i)));
        if(m>0)
            startendlist(:,end+1)=[newend+1;startendlist(2,i);newend+imax];
        end
    end
    
    startendlist(2,i)=newend-1;
    newdeb=startendlist(1,i)-1+find([-1 s(startendlist(1,i):startendlist(3,i))<0],1,'last');
    
    if (newdeb-startendlist(1,i))>=15*4
        [m,imax]=max(s(startendlist(1,i):newdeb-1));
        if(m>0)
            startendlist(:,end+1)=[startendlist(1,i);newdeb-1;startendlist(1,i)+imax-1]; %#ok<*AGROW>            
        end
    end
    
    startendlist(1,i)=newdeb;
    i=i+1;
end
end

%*************************************************************************
function startendlist=accidentcandidat(s1,s2,seuil)
binsig=[0 s1-s2>seuil 0];
startendlist= [ find(binsig(2:end)&~binsig(1:end-1)) ; find(~binsig(2:end)&binsig(1:end-1))-1 ];

for i=1:size(startendlist,2)
    [~,imax]=max(s1(startendlist(1,i):startendlist(2,i))-s2(startendlist(1,i):startendlist(2,i)));
    startendlist(3,i)=startendlist(1,i)+imax-1;
end

end

%*************************************************************************
function [Y,mp]=medgliss(X,win,coef,decim,X2,p2,c)

Xd=butterfilt(X,240,0,240/2.2/decim,8,1);
Xd=Xd(1:decim:end);
if(exist('X2','var'))
    X2=X2(1:decim:end);
end
coefd=decimate(double(coef),decim);
coefd=coefd.*(coefd>0);
midwin=(length(win)-1)/2;
Yd=zeros(size(Xd));

mintolerated=zeros(1,length(Xd));
maxtolerated=255*ones(1,length(Xd));
for i=1:240/2/decim:length(Xd)-10*240/decim
    twin=i:i+10*240/decim;
    mi=min(Xd(twin));
    ma=max(Xd(twin));
    
    mintolerated(twin(mintolerated(twin)<=mi))=mi;
    maxtolerated(twin(maxtolerated(twin)<=ma))=ma;
end

for i=1:length(Xd)
    if(i-1<midwin)
        mwm=i-1;
        mwp=min(max(mwm, floor((midwin-mwm)/2)), length(Xd)-i);
    elseif(length(Xd)-i<midwin)
        mwp=length(Xd)-i;
        mwm=min(max(mwp, floor((midwin-mwp)/2)), i-1);
    else
        mwm=midwin;
        mwp=midwin;
    end
    points=i-mwm:i+mwp;
    Xpoints=Xd(points);
    coefwin=coefd(points).*win(midwin-mwm+1:midwin+mwp+1).*(Xpoints>=mintolerated(i) & Xpoints<=maxtolerated(i));
    s=sum(coefwin);
    scoef=sum(win(midwin-mwm+1:midwin+mwp+1));
    
    if(exist('X2','var'))
        
        Xpoints=[Xpoints X2(i)];
        coefwin=[coefwin max(0,c*p2(i)*scoef-s)];
        s=sum(coefwin);
    end
    [p,order]=sort(Xpoints);
    Yd(i)=p(find(cumsum(coefwin(order))>=s/2,1,'first'));
    mp(i)=s/scoef;
end

Y=interp(Yd,decim);
Y=Y(1:length(X));

mp=interp(mp,decim);
mp=mp(1:length(X));

end


%*************************************************************************
function [y,x]=enveloppe(x,srate,f0,f1)

fftx=fft(x);
siglen=length(x)/srate;
firstsamp=round(f0*siglen+1);
lastsamp=round(f1*siglen+1);

ffty=zeros(size(x));

if(mod(lastsamp-firstsamp,2)==1)
    ffty([end-(lastsamp-firstsamp-3)/2:end 1:(3+lastsamp-firstsamp)/2])=fftx(firstsamp:lastsamp);
    
else
    ffty([end-(lastsamp-firstsamp-2)/2:end 1:(2+lastsamp-firstsamp)/2])=fftx(firstsamp:lastsamp);
end

if(nargout==2)
    fftx([1:firstsamp-1 lastsamp+1:end-lastsamp+1 end-firstsamp+2:end])=0;
    x=ifft(fftx);
end
y=2*abs(ifft(ffty));
end