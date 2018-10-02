% Method T
% Re-coded by S. Boudet from
% Taylor, G. M., Mires, G. J., Abel, E. W., Tsantis, S., Farrell, T., Chien, P. F., Liu, Y. - "The development and validation of an algorithm for real-time computerised fetal heart rate monitoring in labour" BJOG: an international journal of obstetrics and gynaecology 107(9):1130--1137, sep 2000
% 
% USAGE
%    [baseline,accelerations,decelerations]=aamtaylor(FHR)
%         Taylor's method with its method for acceleration/deceleration detection
%    [baseline,accelerations,decelerations]=aamtaylor(FHR,1)
%         Taylor's method with a standard simple method for acceleration/deceleration detection
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
function [baseline,accelerations,decelerations]=aamtaylor(FHR,simpleacc)


FHRi=FHR;
FHRbl=FHRi;

bl1=butterfilt(FHRi,4,0,0.008,3,1);
FHRbl(FHRi-bl1>5)=0;
FHRbl(FHRi-bl1<-5)=0;

bl2=butterfilt(interpolFHR(FHRbl),4,0,0.006,3,1);
FHRbl(FHRi-bl2>5)=0;
FHRbl(FHRi-bl2<-5)=0;

bl3=butterfilt(interpolFHR(FHRbl),4,0,0.006,3,1);
FHRbl(FHRi-bl3>10)=0;
FHRbl(FHRi-bl3<-5)=0;

baseline=butterfilt(interpolFHR(FHRbl),4,0,0.006,3,1);

FuzzyLine=butterfilt(FHRi,4,0,0.02,4,1);

% The following is not clearly described on the paper. 
% The begining and ending of accident is described for the baseline caluclation but
% not directly for the accident determination which we guess should be
% different. Here is the sentence which let us suppose our code.
% 
% "Accelerations and decelerations were determined and classified according
% to definitions given by FIGO" 
% "Artefact lust be removed prior to determining the timing of the
% decelerations; this is achieved by a low pass fourth order Butterworth
% filter with a cut off frequency of 0.02Hz"
% Context on those sentences have ambiguities so there is some doubts

if( exist('simpleacc','var') && simpleacc==1)
    [accelerations,decelerations]=simpleaddetection(FHR,baseline);

else
    %List of periods where there is 10bpm difference which are acceleration
    %candidats
    accelerations=startendlist(FuzzyLine-baseline>10)/4; 
    %Delete periods of too short duration or too low amplitude
    accelerations=validaccident(accelerations,FHRi-baseline,15,15); 


    decelerations=startendlist(FuzzyLine-baseline<-5)/4;
    decelerations=validaccident(decelerations,baseline-FHRi,15,15);
end




