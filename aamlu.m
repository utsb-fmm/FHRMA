% Method L
% Re-coded by S. Boudet from
% LU, Yaosheng et WEI, Shouyi. 
% Nonlinear baseline estimation of FHR signal using empirical mode decomposition. 
% In : Signal Processing (ICSP), 2012 IEEE 11th International Conference on. IEEE, 2012. p. 1645-1649.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aamlu(FHR)
%         Lu's method standard simple method for acceleration/deceleration detection 
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
function  [ baseline ,acc,dec,falseacc,falsedec] =aamlu(FHR0)
j=1;
%Sampling rate ?
l=length(FHR0);
FHR=subsamp(FHR0,4);

d=FHR;
r=FHR;
for j=1:8 %pour l'essai
    SD=1;
    while SD>0.1
        [eu,el]=enveloppe(d);%OK
        m=(eu+el)/2;%OK
        dprev=d;%ok
        d=d-m;%ok
        SD=sum((dprev-d).^2)/sum(dprev.^2);
        %SD=sum((dprev(30:end-30)-d(30:end-30)).^2)/sum(dprev(30:end-30).^2);%30???
    end
    c(j,:)=d;
    r(j+1,:)=r(j,:)-c(j,:); %On the paper it's FHR-c(j,:) but it does not correspond to eq. 2.
    %figure;plot(c(j,:))
    d=FHR-sum(c);
end
baseline=sum(c(3:end,:))+r(j+1,:); 
% this is strange:  baseline can simply be written =FHR-c1-c2

baseline=exclude(baseline,FHR,4);
baseline=resamp(baseline,4,l);

    
[acc,dec,falseacc,falsedec]=simpleaddetection(FHR0,baseline);

end

function sFHR=subsamp(FHR,factor)
    sFHR=zeros(1,floor(length(FHR)/factor));
    for i=1:length(sFHR)
        sFHR(i)=mean(FHR((i-1)*factor+1:i*factor));
    end
end

function [eu,el]=enveloppe(x)
upoints=[1 find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>x(3:end))+1 length(x)];
lpoints=[1 find(x(2:end-1)<=x(1:end-2) & x(2:end-1)<x(3:end))+1 length(x)];
eu = spline(upoints,x(upoints),1:length(x));
el = spline(lpoints,x(lpoints),1:length(x));
end

function baseline=exclude(x,FHR,factor)
    upoints=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>x(3:end))+1 ;
    lpoints=find(x(2:end-1)<=x(1:end-2) & x(2:end-1)<x(3:end))+1 ;
    points=sort([1 upoints lpoints length(x)]);
    D=abs(x(points(2:end))-x(points(1:end-1)));%Simpler but not exactly what is writted
    %Dlim=mean(D(2:end-2))+.5*std(D(2:end-2));
    Dlim=mean(D)+.5*std(D);
    truepoints=D<Dlim;    
    truepoints=[1 truepoints] & [truepoints 1]; % Test with a "or"
    
    CP=[];
    %sliding window of 20min
    winstartlist=0:2*240/factor:length(x)-20*240/factor;
    workinginterval=[winstartlist+9*240/factor ; winstartlist+11*240/factor];
    workinginterval(1,1)=0;workinginterval(2,end)=length(x);
    
    for i=1:length(winstartlist);
        winstart=winstartlist(i);
        p=points(truepoints & points>winstart & points<=winstart+20*240/factor);
        ME=mean(x(p));
        ST=std(x(p));
        pinterval=points(truepoints & points>workinginterval(1,i) & points<=workinginterval(2,i));
        CP=[CP pinterval(abs(x(pinterval)-ME)<ST) ];
        
        
    end
    
%     CP=sort(CP);
%     CP=[CP(CP(1:end-1)<CP(2:end)-1) CP(end)];
    
    % The paper recommend the using of spline which would be
    %baseline=spline(CP,x(CP),1:CP(end));
    %baseline(end:length(x))=baseline(end);
    
    %But the method is very stabler with linear interpolation :
    baseline=linearinterpolation(CP,x(CP),1:length(x));
    %figure;plot((1:length(FHR))/240*factor,[FHR;x;baseline]',points(truepoints)/240*factor,x(points(truepoints)),'*',CP/240*factor,x(CP),'+');
    
    baseline=avgwin(baseline,240/factor); %Win size not given
    
    %me=(x(points(truepoints(2:end)))+x(points(truepoints(1:end-1))))*.5;
    %melim=[mean(le)+std(le) mean(le)-std(le)];
    
    
end

function yy=linearinterpolation(x,y,xx)
    yy=zeros(size(xx));
    yy(1:x(1))=y(1);
    for i=1:length(x)-1
        yy(x(i):x(i+1))=linspace(y(i),y(i+1),x(i+1)-x(i)+1);
    end
    yy(x(end):end)=y(end);
end


