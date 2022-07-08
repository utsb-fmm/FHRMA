% This is the main function to explore the database with the interface
% to annotate FHR Doppler&MHR signals and to display results of FSDop and FSMHR. 
% This is the function to launch if you want to
% explore the database.
% To display a specific recording providing the filename, use
% showFSAnalysis instead. 
% USAGE
%     FSDopGUI
%     FSDopGUI(n)
% INPUT 
%     n    : the number of recording (default 1)


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
classdef FSDopGUI
    properties
        ploter
        listfiles
        model
    end
    methods
        function obj=FSDopGUI(n)
            load('./FSdataset/expertAnnotations.mat','DB');
            DB=DB(contains({DB.dataset},'Dop'));
            obj.listfiles={DB.filename};
            if nargin==0
                [fname,~]=uigetfile('FSdataset/DopMHR/*.fhrm;*.dat');
                n=find(strcmp(fname,obj.listfiles));
            elseif isstring(n)
                [~,fname,ext]=fileparts(n);
                fname=[fname ext];
                n=find(strcmp(fname,obj.listfiles));
            else
                fname=obj.listfiles{n};
            end
            
            d=dir(['FSDataset/**/' fname]);
            if isempty(d)
                errordlg('File must be in FSdataset folder');
                exit()
            end

            obj.ploter=showFSAnalysis([d.folder '/' d.name],'Dop');
            set(obj.ploter.EdtRec,'String',n)
            addlistener(obj.ploter,'ChangeRec',@(src,evt) ChangeRec(obj));
        end

        function fold=DS2fold(~,DS,file)
            if strcmp(DS(7:end),'TestDbS')
                TVT='TestDoubleSignals';
            elseif strcmp(DS(7:end),'TestCP')
                TVT='TestCurrentPractice';
            else
                TVT=DS(7:end);
            end
            fold=[DS(1:6) '/' TVT '/' file];
        end
        
        function ChangeRec(obj)
            n=str2double(get(obj.ploter.EdtRec,'String'));
            close(obj.ploter.Fig)
            delete(obj.ploter);
            d=dir(['FSDataset/**/' obj.listfiles{n}]);
            obj.ploter=showFSAnalysis([d.folder '/' d.name],'Dop');
            set(obj.ploter.EdtRec,'String',n)
            addlistener(obj.ploter,'ChangeRec',@(src,evt) ChangeRec(obj));
        end
    end
end
