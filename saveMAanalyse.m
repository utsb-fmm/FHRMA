%Script which save all method results of the database in a specific file
%for each method

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

load FHRMAdataset/analyses/expertAnalyses.mat
for i=1:156
    if i<=90
        f='FHRMAdataset/testdata/';
    else
        f='FHRMAdataset/traindata/';
    end
    
    
    [FHR1,FHR2,~,TOCO,timestamp]=fhropen([f data(i).filename]);
    [FHRi,FHRraw,TOCOi]=preprocess(FHR1,FHR2,TOCO,data(i).unreliableSignal);
    FHR0=FHRraw;FHR0(isnan(FHR0))=0;
    
     [ WMFB_orig.data(i).baseline,ac,dc ] = aamwmfb(FHRi);
     WMFB_orig.data(i).accelerations=ac(1:2,:)'/60;
     WMFB_orig.data(i).decelerations=dc(1:2,:)'/60;    
     WMFB_orig.data(i).filename=data(i).filename;
    
    
    disp(i)
end

save('myMACustomMethod.mat','-struct','WMFB_orig');

load FHRMAdataset/analyses/expertAnalyses.mat
fileresults2=evaluateaam('myMACustomMethod.mat',data([data(:).trainingData]==1),'FHRMAdataset/traindata/');

