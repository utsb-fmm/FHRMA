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

for i=1:156
    if i<=66
        f='traindata/';
    else
        f='testdata/';
    end
    
    
    [FHR1,FHR2,TOCO,timestamp]=fhropen([f data(i).filename]);
    [~,FHR,TOCOi,d,f]=preprocess(FHR1,FHR2,TOCO);
    rjct=data(i).unreliableSignal;
    for j=1:size(rjct,1)
        FHR(round(rjct(j,1)*240+1):round(rjct(j,2)*240))=NaN;
    end
    FHR0=FHR;FHR0(isnan(FHR0))=0;
    [FHRi,FHR]=preprocess(FHR0,zeros(1,length(FHR0)),TOCO);
    FHR0=FHR;FHR0(isnan(FHR0))=0;
    
    [ T_orig.data(i).baseline,ac,dc ] = aamtaylor(FHRi);
    T_orig.data(i).accelerations=ac(1:2,:)'/60;
    T_orig.data(i).decelerations=dc(1:2,:)'/60;
    
    [ T_std.data(i).baseline,ac,dc] = aamtaylor(FHRi,1);
    T_std.data(i).accelerations=ac(1:2,:)'/60;
    T_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ P_std.data(i).baseline,ac,dc] = aampardey( FHRi );
    P_std.data(i).accelerations=ac(1:2,:)'/60;
    P_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ J_orig.data(i).baseline,ac,dc ] = aamjimenez( FHRi );
    J_orig.data(i).accelerations=ac(1:2,:)'/60;
    J_orig.data(i).decelerations=dc(1:2,:)'/60;
    
    [ J_std.data(i).baseline,ac,dc ] = aamjimenez( FHRi,1 );
    J_std.data(i).accelerations=ac(1:2,:)'/60;
    J_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ H_std.data(i).baseline,ac,dc ] = aamhouze( FHRi,TOCOi ,1);
    H_std.data(i).accelerations=ac(1:2,:)'/60;
    H_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ MD_std.data(i).baseline,ac,dc ] = aammaeda ( FHRi );
    MD_std.data(i).accelerations=ac(1:2,:)'/60;
    MD_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ A_std.data(i).baseline,ac,dc ] = aamayres ( FHRi );
    A_std.data(i).accelerations=ac(1:2,:)'/60;
    A_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ MT_orig.data(i).baseline,ac,dc ] = aammantel ( FHRi ,FHR0);
    MT_orig.data(i).accelerations=ac(1:2,:)'/60;
    MT_orig.data(i).decelerations=dc(1:2,:)'/60;
    
    [ MT_std.data(i).baseline,ac,dc ] = aammantel ( FHRi ,FHR0,1);
    MT_std.data(i).accelerations=ac(1:2,:)'/60;
    MT_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ C_std.data(i).baseline,ac,dc ] = aamcazares ( FHRi );
    C_std.data(i).accelerations=ac(1:2,:)'/60;
    C_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ C_orig.data(i).baseline,ac,dc ] = aamcazares ( FHRi,1 );
    C_orig.data(i).accelerations=ac(1:2,:)'/60;
    C_orig.data(i).decelerations=dc(1:2,:)'/60;
    
    [ L_std.data(i).baseline,ac,dc ] = aamlu ( FHRi );
    L_std.data(i).accelerations=ac(1:2,:)'/60;
    L_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ W_std.data(i).baseline,ac,dc ] = aamwrobel(FHRi);
    W_std.data(i).accelerations=ac(1:2,:)'/60;
    W_std.data(i).decelerations=dc(1:2,:)'/60;
    
    [ MG_std.data(i).baseline,ac,dc ] = aammongelli(FHRi);
    MG_std.data(i).accelerations=ac(1:2,:)'/60;
    MG_std.data(i).decelerations=dc(1:2,:)'/60;
    
    J_orig.data(i).filename=data(i).filename;
    J_std.data(i).filename=data(i).filename;
    H_std.data(i).filename=data(i).filename;
    MD_std.data(i).filename=data(i).filename;
    T_orig.data(i).filename=data(i).filename;
    T_std.data(i).filename=data(i).filename;
    P_std.data(i).filename=data(i).filename;
    MT_orig.data(i).filename=data(i).filename;
    MT_std.data(i).filename=data(i).filename;
    A_std.data(i).filename=data(i).filename;
    C_std.data(i).filename=data(i).filename;
    C_orig.data(i).filename=data(i).filename;
    L_std.data(i).filename=data(i).filename;
    W_std.data(i).filename=data(i).filename;
    MG_std.data(i).filename=data(i).filename;
   
    
    
    
    
    disp(i)
end

save('J_orig','-struct','J_orig');
save('J_std','-struct','J_std');
save('H_std','-struct','H_std');
save('MD_std','-struct','MD_std')
save('T_orig','-struct','T_orig');
save('T_std','-struct','T_std');
save('P_std','-struct','P_std');
save('MT_orig','-struct','MT_orig');
save('MT_std','-struct','MT_std');
save('A_std','-struct','A_std');
save('C_std','-struct','C_std');
save('C_orig','-struct','C_orig');
save('L_std','-struct','L_std');
save('W_std','-struct','W_std');
save('MG_std','-struct','MG_std');