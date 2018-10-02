% Linear Interpollation and extrapollation 
%
% USAGE (same as spline)
%     yy=linearinterpolation(x,y,xx)
%
% INPUT
%     x      : x coordiantes of samples to interpole
%     y      : y coordiantes of samples to interpole
%     xx     : new x coordinates which are get by interpollation
% 
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
function yy=linearinterpolation(x,y,xx)
    yy=zeros(size(xx));
    a=(y(2)-y(1))/(x(2)-x(1));
    b=y(1)-a*x(1);
    yy(xx<x(1))=a*xx(xx<x(1))+b;
    for i=1:length(x)-1
        a=(y(i+1)-y(i))/(x(i+1)-x(i));
        b=y(i)-a*x(i);
        yy(xx>=x(i) &  xx<x(i+1))=a*xx(xx>=x(i) &  xx<x(i+1))+b;
    end
    yy(xx>=x(end))=a*xx(xx>=x(end))+b;
end