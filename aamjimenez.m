% Method J
% Re-coded by S. Boudet from
% JIMENEZ, L., GONZALEZ, R., GAITAN, M., et al. 
% Computerized algorithm for baseline estimation of fetal heart rate. 
% In : Computers in Cardiology, 2002. IEEE, 2002. p. 477-480.
%
% USAGE
%    [baseline,accelerations,decelerations]=aamjimenez(FHR)
%         Jimenez's method with its method for acceleration/deceleration detection
%    [baseline,accelerations,decelerations]=aamjimenez(FHR,1)
%         Jimenez's method with a standard simple method for acceleration/deceleration detection
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


function [baseline,accelerations,decelerations,accelerationsb,decelerationsb]=aamjimenez(FHR,simpleacc)
    FHR0=FHR;
    l=length(FHR);
    FHR=avgsubsamp(FHR,2);

    % Normally sampled as one per beat with thus a sampling rate varying
    % which would make the code complicated
    % To simplify, we consider a sampling rate as 120 samples / min (=240/2)
    % Aplying same parameters as described on the paper will lead to a less
    % smoothed signal on high frequencies.
    srate=2; %Sampling rate in Hz

    %Step 1 smoothing
    sFHR=[FHR(13:-1:1) FHR FHR(end:-1:end-12)];
    h=hann(27);
    sFHR=conv(h/sum(h),sFHR);
    sFHR=sFHR(27:end-26);
    dFHR=(sFHR(2:end)-sFHR(1:end-1))*srate;

    blpoints=abs(dFHR)<=1;
    blpoints(end+1)=blpoints(end);
    for wind=[1:5*60*srate:length(FHR)-5*60*srate length(FHR)-5*60*srate+1 ]
        %Jimenez worked on signal of length 5 min so we decided to divide
        %the signal on 5 min window
        win=wind:wind+5*60*srate-1;
        for i=[1 find(blpoints(win(1:end-1))==0 & blpoints(win(2:end)))]
            %For each begining of stable period
            d=find([blpoints(win(i+1:end)) 0]==0,1,'first');
            if d<15*srate
                blpoints(win(i:i+d-1))=0; %reject to short periods
            end
        end

        % On the original paper :
        % MUwin=mean( sFHR(win( blpoints(win) )) );
        % but problems appear if there is 2 far away periods.
        % The average is then in middle and there is no point
        % in Mu+-10. So median is better appropriated
        MUwin=median( sFHR(win( blpoints(win) )) );

        for i=[1 find(blpoints(win(1:end-1))==0 & blpoints(win(2:end)))]
            d=find([blpoints(win(i+1:end)) 0]==0,1,'first'); %end of stable segment
            MUsegment=mean(sFHR(win(i+1:i+d-1)));
            if abs(MUsegment-MUwin)>10  % Reject far from mean periods
                blpoints(win(i:i+d-1))=0;
            end       
        end
    end

    sFHR(blpoints==0)=0;
    sFHR=interpol(sFHR);

    baseline=butterfilt(sFHR,2,0,0.033,3,1);
    [accelerations,decelerations]=accidents(baseline,FHR);
    accelerationsb=zeros(2,0);
    decelerationsb=zeros(2,0);
    baseline=resamp(baseline,2,l);
    if( exist('simpleacc','var') && simpleacc==1)
        [accelerations,decelerations,accelerationsb,decelerationsb]=simpleaddetection(FHR0,baseline);
    end
end

function A=interpol(B)
    D=find(B(1:end-1)>0 & B(2:end)==0);
    F=find(B(1:end-1)==0 & B(2:end)>0);
    A=B;
    A(B==0)=0;
    if(F(1)<D(1))
        delta=(B(F(1)+1)-B(F(1)+11))/10; %The derivative is computed on 5s (completely arbitrary)
        A(1:F(1))=linspace(B(F(1)+1)+F(1)*delta,B(F(1)+1)+delta,F(1));
        F=F(2:end);
    end

    if(D(end)>F(end))
        delta=(B(D(end))-B(D(end)-10))/10;
        A(D(end)+1:end)=linspace(B(D(end))+delta,B(D(end))+(length(A)-D(end))*delta,length(A)-D(end));
        D=D(1:end-1);
    end

    for i=1:length(D)
        x=[D(i)-10 D(i) F(i)+1 F(i)+10];
        y=B(x);
        A(D(i)+1:F(i)) = spline(x,y,D(i)+1:F(i));
    end
end

function [accelerations,decelerations]=accidents(BL,FHR)
    diff=BL-FHR;
    transpoint=[1 find(diff(1:end-1)>=0 & diff(2:end)<0 | diff(1:end-1)<0 & diff(2:end)>=0) length(FHR)];
    accelerations=zeros(2,0);
    decelerations=zeros(2,0);
    for i=1:length(transpoint)-1
        t=transpoint(i)+1:transpoint(i+1);
        if(diff(transpoint(i)+1)>=0)
            [m,p]=max(diff(t));
            if length(t)>=30 && p<60 && m>=15
                decelerations=[decelerations [transpoint(i)+1 ;transpoint(i+1)]/2];
            end
        else
            [m,p]=min(diff(t));
            if length(t)>=30 && p<60 && m<=-15
                accelerations=[accelerations [transpoint(i)+1 ;transpoint(i+1)]/2];
            end
        end
    end
end

