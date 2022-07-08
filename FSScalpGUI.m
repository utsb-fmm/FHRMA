% This is the main function to explore the database with the interface
% to annotate FHR Scalp ECG signals and to display results of FSScalp. 
% This is the function to launch if you want to
% explore the database.
% To display a specific recording providing the filename, use
% showFSAnalysis instead. 
% USAGE
%     FSScalpGUI
%     FSScalpGUI(n)
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
classdef FSScalpGUI
    properties
        ploter
        listfiles
        model
    end
    methods
        function obj=FSScalpGUI(n)
            load('./FSdataset/expertAnnotations.mat','DB');
            DB=DB(contains({DB.dataset},'Dop'));
            obj.listfiles={DB.filename};
            if nargin==0
                [fname,~]=uigetfile('FSdataset/ScalpECG/*.fhrm;*.dat');
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

            obj.ploter=showFSAnalysis([d.folder '/' d.name],'Scalp');
            set(obj.ploter.EdtRec,'String',n)
            addlistener(obj.ploter,'ChangeRec',@(src,evt) ChangeRec(obj));
        end

        function fold=DS2fold(~,DS,file)
            fold=['ScalpECG/' DS(6:end) '/' file];
        end
        
        function ChangeRec(obj)
            n=str2double(get(obj.ploter.EdtRec,'String'));
            close(obj.ploter.Fig)
            delete(obj.ploter);
            d=dir(['FSDataset/**/' obj.listfiles{n}]);
            obj.ploter=showFSAnalysis([d.folder '/' d.name],'Scalp');            
            set(obj.ploter.EdtRec,'String',n)
            addlistener(obj.ploter,'ChangeRec',@(src,evt) ChangeRec(obj));
        end
    end
end