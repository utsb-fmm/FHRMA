% BUTTERFILT Complete butterworth filter process 
%   USAGE :
%  >>data=butterfilt(data,srate,f1,f2,order,zeroPhase)
%
% INPUT:
%  data      :the multi channel signals. Signals have to be on lines
%  srate     :Sampling rate
%  f1, f2,   :Frequency cutoff bands
%             if f1=0   Low pass at f2
%             if f2=0   High pass at f1
%             if f1<f2  Band Pass at [f1 f2]
%             if f1>f2  Band stop at [f2 f1]
%  Order     :Order of the filter
%  zeroPhase :<true|false> true = Forward-Backward Filter

%123456789012345678901234567890123456789012345678901234567890123456789012
%
%     BioSigPlot Copyright (C) 2010 Samuel Boudet, Faculté Libre de Médecine,
%     samuel.boudet@gmail.com
%
%     This file is part of BioSigPlot
%
%     BioSigPlot is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     BioSigPlot is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% V0.1 Beta - 24/03/2010 - Initial Version

function data=butterfilt(data,srate,f1,f2,order,zeroPhase)
if ~exist('order','var'), order=6; end
if ~exist('zeroPhase','var'), zeroPhase=1; end

if f1>0 && f2==0
    [b1,a1] = butter(order,2*f1/srate,'high');
elseif f1==0 && f2>0
    [b1,a1] = butter(order,2*f2/srate,'low');
elseif f1>0 && f2>0 && f2<f1
    [b1,a1] = butter(order,2*[f2 f1]/srate,'stop');
elseif f1>0 && f2>0 && f2>f1
    [b1,a1] = butter(order,2*[f1 f2]/srate);
end

if zeroPhase
    nfilt = max(length(b1),length(a1));
    
    if length(b1) < nfilt, b1(nfilt)=0; end   % zero-pad if necessary
    if length(a1) < nfilt, a1(nfilt)=0; end

    rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
    cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
    tdata = [1+a1(2) a1(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
    sp = sparse(rows,cols,tdata);
    zi = sp \ ( b1(2:nfilt).' - a1(2:nfilt).'*b1(1) );
    
    if strcmp(class(data),'single')
        data=single(multisigfilter(b1,a1,double(data'),zi))';
    else
        data=multisigfilter(b1,a1,data',zi)';
    end
else
    
    for i=1:size(data,1)
        if strcmp(class(data),'single')
            data(i,:)=single(filter(b1,a1,double(data(i,:))));
        else
            data(i,:)=filter(b1,a1,data(i,:));
        end
    end
end

