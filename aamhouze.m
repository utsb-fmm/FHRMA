% Method H
% Re-coded by S. Boudet from original program published in
% D. Houzé de l’Aulnoit, R. Beuscart, G. Brabant, L. Carette, M. Delcroix, 
% Real-time Analysis Of The Fetal Heart Rate, 
% in: IEEE, 1990: pp. 1994–1995.
%
%
% USAGE
%    [baseline,accelerations,decelerations]=aamhouze(FHR,TOCO,1)
%         Houzé's method with standard simple method for acceleration/deceleration detection
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
%     FHR Morphological Analysis Toolbox Copyright (C) 2018 Samuel Boudet, Faculté de Médecine et Maïeutique,
%     samuel.boudet@gmail.com
%
%     This file is part of FHR Morphological Analysis Toolbox
%
%     FHR Morphological Analysis Toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     FHR Morphological Analysis Toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [baseline,accelerations,decelerations,accelerationsb,decelerationsb]=aamhouze(FHR,TOCO,simpleacc)

[BL,AccelInfor,RalenInfor] = FHRanalysor(FHR, TOCO);

baseline=zeros(size(FHR));

for i=1:size(BL,2)-1
    baseline(BL(4,i):BL(4,i+1))=linspace(BL(2,i),BL(2,i+1),BL(4,i+1)-BL(4,i)+1);
end
baseline(1:BL(4,3))=BL(2,3);
baseline(BL(4,end):end)=BL(2,end);

if( exist('simpleacc','var') && simpleacc==1)
    [accelerations,decelerations,accelerationsb,decelerationsb]=simpleaddetection(FHR,baseline);
else
    
    accelerations=reshape(AccelInfor(2:floor((end-1)/2)*2+1),2,floor((length(AccelInfor)-1)/2))/4;
    decelerations=reshape(RalenInfor{1}(1:floor(end/2)*2),2,floor(length(RalenInfor{1})/2))/4;
    

end

end

function [FHR_BaseLine,AccelInfor,RalenInfor] = FHRanalysor(FHR1, TOCO)

FHR_Lys = CourtLys(FHR1);% lissage court
[CuInfor, CuState, TOCO_Lys] = DetectCu(TOCO);
% CuInfor: 1*3 matrice, 1 for mark, 2 for detail ampli, 3 for top
% CuState: utile dans la fonction BaseLine, pour determiner un ralentissement
% TOCO_Lys: lissage court
[FHR_LongLys, diff_coef] = MeanRythme(FHR1);
% FHR_LongLys: utile dans la fonciton baseline
% diff_coef: 1*2 matrice, 1 pour derive de rcflis, 2 pour derive secondaire de rcflis
[AccelInfor, RalenInfor, FHR_BaseLine] = BaseLine(FHR_Lys, CuState, FHR_LongLys, diff_coef);
% AccelInfor: marque des accelerations
% RalenInfor: 1*5 matrice, 1 for mark, 2 for type, 3 for detail, 4 for color, 5 for danger
% FHR_BaseLine: 1*4 matrice, 1 for sup, 2 for mean, 3 for inf, 4 for time

end

function [infor, CuState, Lys] = DetectCu(input)
% Infor: 1*3 matrice, 1 for mark, 2 for detail ampli, 3 for top
% CuState: utile dans la fonction BaseLine, pour determiner un ralentissement
% Lys: lissage court

p = mean(reshape(input(1:180),12,length(input(1:180))/12));

for i = 1 : 11
    l(i) = (p(i) + p(i+1) + p(i+3) + p(i+4))/4;
    p(i+2) = l(i);
end
a = [-21 5 26 41 50 53 50 41 26 5 -21]*l' /256;
if abs(l(6) - a) > a*0.05
    l(6) = a;
end
Lys(1) = l(6);

% structure cu
cu.diff_cu = 0;
cu.cont = 0;
cu.ampli_Max_cu = 0;
cu.Min_d = 0;
cu.fin = 0;
cu.fincu = 1;
cu.debutcu = 0;
cu.seuil = 0;
cu.Max_cu = 0;
cu.count = [0 0 0];
infor = 0;
CuState = [];

for PROCESS_count = 192 : 12 : length(input)
    l(1 : 10) = l(2 : 11);
    p = mean(reshape(input(PROCESS_count-23 : PROCESS_count),12,2));
    l(11) = (l(9) + l(10) + p(1) + p(2))/4;
    a = [-21 5 26 41 50 53 50 41 26 5 -21]*l' /256;
    if abs(l(6) - a) > a*0.05
        l(6) = a;
    end
    Lys(floor((PROCESS_count-168)/12)) = l(6);
    
    % enregistrement de CuState
    if (PROCESS_count - 180)/12 < 65
        
        if cu.fincu == 0
            if cu.debutcu > 0 && cu.debutcu < 6
                CuState(floor((PROCESS_count - 180)/12)) = 1;
            end
            if cu.debutcu > 5
                CuState(floor((PROCESS_count - 180)/12)) = 2;
            end
        else
            CuState(floor((PROCESS_count - 180)/12)) = 3;
        end
        
    else
        if cu.fincu == 0
            if cu.debutcu > 0 && cu.debutcu < 11
                CuState(floor((PROCESS_count - 180)/12)) = 1;
            end% 5*4 * .8%%=16%%/debut.cu=precoce
            if cu.debutcu > 10
                CuState(floor((PROCESS_count - 180)/12)) = 2;
            end
        else
            CuState(floor((PROCESS_count - 180)/12)) = 3;
        end
    end
    
    [cu infor] = Contraction(cu, infor, PROCESS_count, l);% detection de contraction
    
end


end


function [rcflis, diff_coef] = MeanRythme(input)
% rcflis: utile dans la fonciton baseline
% diff_coef: 1*2 matrice, 1 pour derive de rcflis, 2 pour derive secondaire de rcflis

lc = mean(reshape(input(1:528),48,length(input(1:528))/48));
ac = [-21 5 26 41 50 53 50 41 26 5 -21] * lc' / 256;
diff_coef = {[] []};

if abs(lc(6) - ac) > ac*0.05
    lc(6) = ac;
end
rcflis(1) = lc(6);
diff_coef{1}(1) = [27 18 9]*(lc(9:-1:7) - lc(3:5))'/256;
diff_coef{2}(1) = 0;

for i = 576 : 48 : length(input)
    lc(1:10) = lc(2:11);
    lc(11) = mean(input(i-47 : i));
    ac = [-21 5 26 41 50 53 50 41 26 5 -21] * lc' / 256;
    diff_coef{1}(floor((i-480)/48)) = [27 18 9]*(lc(9:-1:7) - lc(3:5))'/256;
    %     bc = (-25.3125 * lc(3) - 16.875 * lc(4) - 8.4375 * lc(5) + 8.4375 * lc(7) + 16.875 * lc(8) + 25.3125 * lc(9))/256;
    diff_coef{2}(floor((i-480)/48))= diff_coef{1}(floor((i-480)/48)) - diff_coef{1}(floor((i-528)/48));
    if abs(lc(6) - ac) > ac*0.05
        lc(6) = ac;
    end
    rcflis(floor((i-480)/48)) = lc(6);
    %     [DUREE_moyen,DUREE_moyen2] = DUREE_moyenne(rcflis, DUREE_moyen,
    %     DUREE_moyen2);;
end
end

function [accel_infor, ralen_infor, baseLine] = BaseLine(FHR_Lys, CuState, FHR_LongLys, diff_coef)
% accel_infor: marque des accelerations
% ralen_infor: 1*5 matrice, 1 for mark, 2 for type, 3 for detail, 4 for color, 5 for danger
% baseLine: 1*4 matrice, 1 for sup, 2 for mean, 3 for inf, 4 for time

accel_infor = 0;
ralen_infor = {[] [] [] {} []};% 1 for mark, 2 for type, 3 for detail, 4 for color, 5 for danger
baseLine = 0;% 1 for sup, 2 for mean, 3 for inf, 4 for time

% parametre baseline
bl.ectlissupDETECT = 1.1; %1.1 cf. ZAPLACC
bl.ectlisinfDETECT = 1.25; %1.2
bl.slismoyenne = 0;
bl.ect_lis = nan;
bl.ect_lis_sup = 200;
bl.ect_lis_inf = 60;
bl.t_lis_moy = 480;
bl.bclis = zeros(1,11);
bl.bc3lis = zeros(1,11);
bl.slislong = zeros(1,11);
bl.fenetre_rcf = zeros(1,65);
bl.fenetre_cu = zeros(1,65);
% 1 for accel, 2 for ralenMark, 3 for ralenDetail, 4 for danger
bl.count = [2 1 0 0 0];

% Parametre ralentissement
bl.GRILLE = nan;
bl.couleur = [1 1 1];
bl.diff_rcf = 0;
bl.t_tal = 0;
bl.apres_RALEN = nan;
% bl.apres_accel = nan;
bl.DUREE_zone_residuelle = nan;
bl.DUREE_RALENV = 60;
bl.surf_ralenV = 50;
bl.RALEN = 0;
bl.DUREE_RALEN = 0;
bl.debut_RALEN_vari = 0;
bl.RALEN_cu = 0;
bl.surf_RALEN = 0;
bl.bool_RALEN = 0;
bl.Min_RALEN = 0;
bl.aprescu = 0;

% Parametre acceleration
bl.DUREE_accel = 0;
bl.bool_accel = 0;
bl.Max_accel = 0;

% Adaptation des codes en basic
FHR_Lys(67 : end + 1) = FHR_Lys(66 :  end);
FHR_Lys(66) = 0;

for k = 11 : length(FHR_LongLys)   %begin from the 1008th point
    
    bl.slislong = FHR_LongLys(k-10 : k);
    bl.bclis = diff_coef{1}(k-10 : k);
    bl.bc3lis = diff_coef{2}(k-10 : k);
    bl.fenetre_rcf = FHR_Lys(4*k - 38 : 4 * k + 26);
    %the first CuState:£¨1008 - 12 - 180£© / 12 = 68
    bl.fenetre_cu = CuState(4*k - 40 : 4 * k + 24);
    bl.t_lis_moy = bl.t_lis_moy + 48;
    %---------  ralentissement  -----------%
    if bl.RALEN == 1
        bl.apres_RALEN = 0;
        if bl.ect_lis * bl.ectlisinfDETECT < 4
            DUREE_ralenDETECT = 360;
        else
            DUREE_ralenDETECT = 240;
        end
        if (bl.RALEN_cu ~= 1 || bl.DUREE_zone_residuelle > 240) && bl.DUREE_RALEN > DUREE_ralenDETECT && bl.bclis(6) < 1 && bl.bc3lis(6) < 0.5
            if bl.slislong(6) > 100
                if bl.surf_RALEN > 300 && (bl.fenetre_rcf(65) > (bl.ect_lis_inf - 20) || bl.Min_RALEN > 40)
                    % validation
                    bl.count(2) = bl.count(2) + 1;
                    ralen_infor{1}(bl.count(2)) = bl.t_tal;
                    bl.count(2) = bl.count(2) + 1;
                    
                    if bl.debut_RALEN_vari == 0
                        switch bl.RALEN_cu
                            case 1
                                ralen_infor{2}(floor(bl.count(2)/2)) = 1;
                            case 2
                                ralen_infor{2}(floor(bl.count(2)/2)) = 2;
                        end
                    else
                        ralen_infor{2}(floor(bl.count(2)/2)) = 3;
                    end
                end
                %--------------------------------------------------
                % GOSUB RALEN_MAZ;
                %--------------------------------------------------
                [bl] = RALEN_MAZ(bl);
                bl.slislong(1:6) = bl.slislong(7);
                %------------------------------------------------------
                %    GoTo CALCUL;
                %--------------------------------------------------
                if bl.bool_accel == 1
                    [bl, accel_infor] = ACCELERE(bl, accel_infor);
                end% pour mettre fin d'une acceleration
                bl = RECALCUL(bl);
            end
        end
        [bl, ralen_infor] = RALENTISS(bl, ralen_infor);
    elseif bl.slislong(6) < bl.ect_lis_inf && bl.bclis(6) < -0.5 %ne pas toucher DEMO3
        bl.RALEN = 1;
        [bl, ralen_infor] = RALENTISS(bl, ralen_infor);
    end
    
    %---------  acceleration  -----------%
    if bl.bool_accel == 1
        [bl, accel_infor] = ACCELERE(bl, accel_infor);
    elseif bl.slislong(6) > bl.ect_lis_sup
        bl.bool_accel = 1;
        [bl, accel_infor] = ACCELERE(bl, accel_infor);
        
    end
    
    
    %---------  adaptation mean rythme  -----------%
    % aux fluctuations lentes du rythme de base  cf. PLAT
    if bl.slislong(6) > bl.ect_lis_inf && bl.slislong(6) < bl.ect_lis_sup && bl.bc3lis(7) > -0.3 && bl.bc3lis(8) < 0.4 && bl.bc3lis(8) > -0.3 && bl.bool_RALEN == 0
        if bl.bool_accel == 1
            [bl, accel_infor] = ACCELERE(bl, accel_infor);
        end% pour mettre fin d'une acceleration
        bl = RECALCUL(bl);
    end
    % ?la mont‚e
    if bl.slislong(6) > bl.ect_lis_sup && bl.bclis(6) > 0 && bl.RALEN == 0
        if bl.bool_accel == 1
            [bl, accel_infor] = ACCELERE(bl, accel_infor);
        end% pour mettre fin d'une acceleration
        bl = RECALCUL(bl);
    end
    % recalcul ?la fin des accidents
    if bl.apres_RALEN == 1 && bl.bclis(6) > -1
        bl.apres_RALEN = 0;
        if bl.bool_accel == 1
            [bl, accel_infor] = ACCELERE(bl, accel_infor);
        end% pour mettre fin d'une acceleration
        bl = RECALCUL(bl);
    end
    %     if bl.apres_accel == 1
    %         bl.apres_accel = 0;
    %         if bl.bool_accel == 1
    %             [bl, accel_infor] = ACCELERE(bl, accel_infor);
    %         end% pour mettre fin d'une acceleration
    %         bl = RECALCUL(bl);
    %     end
    %  si aucune condition n'est remplie
    
    %-------  PILE LISSAGE LONG  ------------%
    if bl.slismoyenne > 0
        bl.count(5) = bl.count(5) + 1;
        baseLine(1,bl.count(5)) = bl.ect_lis_sup;
        baseLine(2,bl.count(5)) =  bl.slismoyenne;
        baseLine(3,bl.count(5)) =  bl.ect_lis_inf;
        baseLine(4,bl.count(5)) =  bl.t_lis_moy;
    else
        bl.count(5) = bl.count(5) + 1;
        baseLine(1,bl.count(5)) = 140;
        baseLine(2,bl.count(5)) =  140;
        baseLine(3,bl.count(5)) =  140;
        baseLine(4,bl.count(5)) =  bl.t_lis_moy;
    end
    % on bloquee cartmoyen
    FHR_LongLys(k-10 : k) = bl.slislong;
    
    if k  == length(FHR_LongLys)
        if bl.RALEN == 1
            if (bl.DUREE_RALEN > bl.DUREE_RALENV && bl.surf_RALEN > bl.surf_ralenV)     %conditions de validation d%un RALEN
                % Validation
                bl.count(2) = bl.count(2) + 1;
                ralen_infor{1}(bl.count(2)) = bl.t_tal;
                bl.count(2) = bl.count(2) + 1;
                
                if bl.debut_RALEN_vari == 0
                    switch bl.RALEN_cu
                        case 1
                            ralen_infor{2}(floor(bl.count(2)/2)) = 1;
                        case 2
                            ralen_infor{2}(floor(bl.count(2)/2)) = 2;
                    end
                else
                    ralen_infor{2}(floor(bl.count(2)/2)) = 3;
                end
                
            end              %FIN VALIDATION du RALEN
        end
    end
end
end

function Lys = CourtLys(input)

p = mean(reshape(input(1:180),12,length(input(1:180))/12));
for i = 1 : 11
    l(i) = (p(i) + p(i+1) + p(i+3) + p(i+4))/4;
    p(i+2) = l(i);
end
a = [-21 5 26 41 50 53 50 41 26 5 -21]*l' /256;
if abs(l(6) - a) > a*0.05
    l(6) = a;
end
Lys(1) = l(6);

for PROCESS_count = 192 : 12 : length(input)
    l(1 : 10) = l(2 : 11);
    p = mean(reshape(input(PROCESS_count-23 : PROCESS_count),12,2));
    l(11) = (l(9) + l(10) + p(1) + p(2))/4;
    a = [-21 5 26 41 50 53 50 41 26 5 -21]*l' /256;
    if abs(l(6) - a) > a*0.05
        l(6) = a;
    end
    Lys(floor((PROCESS_count-168)/12)) = l(6);
end
end



%--------------------------------------------------------------------------
%               ---------  fonction Contraction  ----------
%--------------------------------------------------------------------------
function [cu infor] = Contraction(cu, infor, tempsa, l)
b = [27 18 9]*(l(9:-1:7) - l(3:5))'/256;
tempsa = tempsa  - 90;
if cu.cont == 0                  %INTERPRETATION DES CONTRACTIONS
    if b > 1                      %DEBUT DE LA CU
        cu.cont = 1;
        cu.seuil = l(5);
        cu.debutcu = 1; cu.Max_cu = 0;
        cu.count(1) = cu.count(1) +1;
        infor(1, cu.count(1)) = tempsa;
        cu.fincu = 0;
    else
        return
    end
else
    %--   valeurs de B>0
    if cu.fin > 0 && b > 1
        cu.seuil = l(5) - 5; cu.Min_d = 0; cu.fin = 0;
    end
    if b > 0
        cu.debutcu = cu.debutcu + 1;      %periode debut de CU de 7*3,2%%
        cu.diff_cu = l(6) - cu.seuil;
        if cu.diff_cu >= 0
            cu.count(2) = cu.count(2) +1;
            infor(4, cu.count(2)) = tempsa;
            infor(2, cu.count(2)) = cu.diff_cu;
            [cu, infor] = TopCu(cu, infor, b, tempsa);
        end
        return
    end
    
    %--   valeurs de B<0 ; derivee minimale non trouvee
    if b < 0 && b < cu.Min_d
        cu.Min_d = b;
        cu.fin = l(6);
        cu.diff_cu = l(6) - cu.seuil;
        if cu.diff_cu >= 0
            cu.count(2) = cu.count(2) +1;
            infor(4, cu.count(2)) = tempsa;
            infor(2, cu.count(2)) = cu.diff_cu;
            [cu, infor] = TopCu(cu, infor, b, tempsa);
        end
        return
    end
    %--   valeurs de B<0 ; derivee minimale trouvee
    if b < 0 && b > cu.Min_d
        cu.fin = cu.fin + cu.Min_d;
    end
    %--   condition de fin de contraction
    if cu.fin <= cu.seuil + 3
        cu.cont = 0; cu.Min_d = 0; cu.debutcu = 0;
        cu.fin = 0;
        cu.fincu = 1;
        %graphique de la fin de la CU
        cu.count(1) = cu.count(1) +1;
        infor(1, cu.count(1)) = tempsa;
        return
    end
    %--   aucune condition remplie
    cu.diff_cu = cu.fin - cu.seuil;
    if cu.diff_cu >= 0
        cu.count(2) = cu.count(2) +1;
        infor(4, cu.count(2)) = tempsa;
        infor(2, cu.count(2)) = cu.diff_cu;
        [cu, infor] = TopCu(cu, infor, b, tempsa);
    end
    return
end
end


function [cu, infor] = TopCu(cu, infor, b, tempsa)
if cu.diff_cu > cu.Max_cu && b < 0.5
    cu.Max_cu = round(cu.diff_cu);
    cu.ampli_Max_cu = 1;
elseif cu.ampli_Max_cu == 1
    cu.count(3) = cu.count(3) +1;
    infor(3, cu.count(3)) = tempsa - 12;
    cu.ampli_Max_cu = 0;
end
end


function bl = RECALCUL(bl) %------  RECALCUL DU RYTHME MOYEN

bl.slismoyenne = 0; slis2 = 0;
for lis = 1 : 11
    bl.slismoyenne = bl.slismoyenne + bl.slislong(lis); slis2 = slis2 + (bl.slislong(lis) ^ 2);
end
bl.ect_lis = sqrt((slis2 - (bl.slismoyenne ^ 2) / 11) / 10);
bl.slismoyenne = bl.slismoyenne / 11;
bl.ect_lis_sup = bl.slismoyenne + (bl.ect_lis * bl.ectlissupDETECT);
bl.ect_lis_inf = bl.slismoyenne - (bl.ect_lis * bl.ectlisinfDETECT);
end


%--------------------------------------------------------------------------------------------
%                                    ralentissement
%            Couleur: PRECOCE: [0.69 0.88 0.9]; TARDIF: [0 0 0]; VARIABLE: [0 1 0];
%--------------------------------------------------------------------------------------------
function [bl, ralen_infor] = RALENTISS(bl, ralen_infor)
%CALCUL DES PARAMETRES DES RALENTISSEMENTS  calcul sur lissage court
for un = 30 : 33
    if bl.DUREE_RALEN <= 120
        bl.GRILLE = 0;
    else
        bl.GRILLE = 3;
    end
    if bl.fenetre_rcf(un) < (bl.ect_lis_inf - bl.GRILLE) || (bl.fenetre_rcf(un) <= bl.slismoyenne && bl.slislong(6) < (bl.ect_lis_inf - bl.GRILLE))
        bl.bool_RALEN = 1;
        if bl.DUREE_RALEN == 12
            ralen_infor{1}(bl.count(2)) = bl.t_tal;
        end
        bl.DUREE_RALEN = bl.DUREE_RALEN + 12;
        bl.diff_rcf = round(bl.slismoyenne - bl.fenetre_rcf(un));
        bl.Min_RALEN = max(bl.diff_rcf, bl.Min_RALEN);
        %  1008 points de retard: 1008 = 48*11 + 48*10 /RCF//////180+12*33
        bl.surf_RALEN = bl.surf_RALEN + bl.diff_rcf;
        %----- r‚partition des diff_rcf‚rentes surfaces
        if bl.DUREE_RALEN == 12                           % type de ralentissement
            switch bl.fenetre_cu(un)
                case 1                                  %RALENTISSEMENT PENDANT CU , PRECOCE
                    bl.RALEN_cu = 1;
                    bl.couleur = [0.69 0.88 0.9];
                    if bl.diff_rcf > 60
                        bl.count(4) = bl.count(4) + 1;
                        ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                    end
                case 2                                  % RALENTISSEMENT PENDANT CU , TARDIF
                    bl.RALEN_cu = 2;
                    bl.couleur = [0 0 0];
                    if bl.diff_rcf > 60
                        bl.count(4) = bl.count(4) + 1;
                        ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                    end
                case 3                                  %  RALENTISSEMENT VARIABLE
                    bl.aprescu = 1;
                    bl.RALEN_cu = 0;
                    bl.debut_RALEN_vari = 1;
                    bl.couleur = [0 1 0];
                    if bl.diff_rcf > 14
                        bl.count(4) = bl.count(4) + 1;
                        ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                    end
            end
        elseif bl.DUREE_RALEN > 12
            if bl.fenetre_cu(un) == 3
                bl.aprescu = 1;
                if bl.RALEN_cu > 0         %  RALENTISSEMENT PENDANT CU , AVEC ZONE RESIDUELLE
                    bl.DUREE_zone_residuelle = bl.DUREE_zone_residuelle + 12;
                    bl.couleur = [0 0 1];
                    if bl.diff_rcf > 30
                        bl.count(4) = bl.count(4) + 1;
                        ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                    end
                elseif bl.RALEN_cu == 0                            %  RALENTISSEMENT VARIABLE
                    bl.couleur = [0 1 0];
                    if bl.diff_rcf > 14
                        bl.count(4) = bl.count(4) + 1;
                        ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                    end
                end
            elseif bl.fenetre_cu(un) < 3
                if bl.aprescu == 1
                    switch bl.debut_RALEN_vari
                        case 1
                            if bl.DUREE_RALEN < 180
                                bl.RALEN_cu = 1;
                                bl.couleur = [0.69 0.88 0.9];
                                bl.debut_RALEN_vari = 0;
                                bl.aprescu = 0;
                            end
                            if bl.diff_rcf > 14
                                bl.count(4) = bl.count(4) + 1;
                                ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                            end
                        case 0                           % RALENTISSEMENT CU mais INTERPRETE COMME VARIABLE
                            bl.DUREE_zone_residuelle = 0;
                            bl.couleur = [0 1 0];
                            bl.debut_RALEN_vari = 1;
                            bl.RALEN_cu = 0;
                            if bl.diff_rcf > 14
                                bl.count(4) = bl.count(4) + 1;
                                ralen_infor{5}(bl.count(4)) =  bl.t_tal;
                            end
                    end
                end
            end
        end
        %----- enregistrement du ralen
        bl.t_tal = bl.t_lis_moy + 12 * (un -33);
        if bl.ect_lis_inf > bl.fenetre_rcf(un)
            bl.count(3) = bl.count(3) + 1;
            ralen_infor{3}(bl.count(3)) = bl.t_tal;
            ralen_infor{4}{bl.count(3)} = bl.couleur;
        end
    else
        if bl.bool_RALEN == 1
            if (bl.DUREE_RALEN > bl.DUREE_RALENV && bl.surf_RALEN > bl.surf_ralenV)     %conditions de validation d%un RALEN
                % validation
                bl.count(2) = bl.count(2) + 1;
                ralen_infor{1}(bl.count(2)) = bl.t_tal;
                bl.count(2) = bl.count(2) + 1;
                
                if bl.debut_RALEN_vari == 0
                    switch bl.RALEN_cu
                        case 1
                            ralen_infor{2}(floor(bl.count(2)/2)) = 1;
                        case 2
                            ralen_infor{2}(floor(bl.count(2)/2)) = 2;
                    end
                else
                    ralen_infor{2}(floor(bl.count(2)/2)) = 3;
                end
                
            end              %FIN VALIDATION du RALEN
            
            %--------------------------------------------------
            % GOSUB RALEN_MAZ;
            %--------------------------------------------------
            [bl] = RALEN_MAZ(bl);
            bl.DUREE_zone_residuelle = 0;
            bl.apres_RALEN = 1;
            bl.slislong(1:6) = bl.slislong(7);
        end
        return
    end
end
end

function [bl] = RALEN_MAZ(bl)
bl.RALEN = 0;
bl.bool_RALEN = 0;
bl.DUREE_RALEN = 0;
bl.surf_RALEN = 0;
bl.debut_RALEN_vari = 0;
bl.RALEN_cu = 0;
bl.couleur = [1 1 1];
bl.Min_RALEN = 0;
bl.aprescu = 0;
end


%--------------------------------------------------------------------------
% acceleration
%--------------------------------------------------------------------------
function [bl, accel_infor] = ACCELERE(bl, accel_infor)

for un = 30 : 33
    if bl.fenetre_rcf(un) > bl.ect_lis_sup || (bl.fenetre_rcf(un) > bl.slismoyenne && bl.slislong(6) > bl.ect_lis_sup)
        if bl.DUREE_accel == 12
            accel_infor(bl.count(1)) = bl.t_tal; % debut d'une acceleration
        end
        bl.DUREE_accel = bl.DUREE_accel + 12;
        bl.diff_rcf  = bl.fenetre_rcf(un) - bl.slismoyenne;
        bl.Max_accel = max(bl.Max_accel, bl.diff_rcf);
        bl.t_tal = bl.t_lis_moy + 12 * (un -33);
    else
        if bl.slislong(6) < bl.ect_lis_sup
            if bl.Max_accel > 3 * bl.ect_lis && bl.Max_accel > 5
                bl.count(1) = bl.count(1) + 1;
                accel_infor(bl.count(1)) = bl.t_tal; % fin d'une acceleration
                bl.count(1) = bl.count(1) + 1;
                % ---- ACCELERATION NON VALIDEE MAIS VALIDATION OSCILLATION
            elseif bl.Max_accel > 1.5 * bl.ect_lis && bl.Max_accel > 5
                %                     fin_oscillation = 1;
                accel_infor(1) = accel_infor(1) + 1;
                %                     DUREE_varimoy(3) = DUREE_varimoy(3) + 1;
            end
            bl.DUREE_accel = 0;
            bl.bool_accel = 0;
            bl. Max_accel = 0;
            return
        end
    end
end
end

