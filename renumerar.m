%  <CÁLCULO DO FLUXO DE POTÊNCIA - LOAD FLOW APPLYING NEWTON-RAPHSON METHOD V1.0. 
%  This is the main source of this software that calculates the power flow of a power network described using an excel input data file >
%     Copyright (C) <2014>  <Sebastián de Jesús Manrique Machado>   <e-mail:sebajmanrique747@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%Cálculo_Fluxo_de_Potência_ Newton Raphson
%   Sebastián de Jesús Manrique Machado
%   Estudante_Mestrado Em Engenharia Elétrica
%   UEL - 2014.


function [ Dados_linhas, Dados_barras, vetor_aux2, num_barras, num_linhas, num_n_PV, num_n_PQ, num_n_novos ] = renumerar(  Dados_linhas, Dados_barras, num_linhas, num_barras )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%||  Ordenamento ||
%==================

aux_novo_n = max([Dados_linhas(:,1) ; Dados_linhas(:,2)].') + 1;
Index_n_Sl  = find(Dados_barras(:,2)==1); 
Index_n_PV  = find(Dados_barras(:,2)==2);                       % posição Barra PV
Index_n_PQ  = find(Dados_barras(:,2)==3);                       % Posição Barra PQ
num_n_PV    = length(Index_n_PV);                               % Quantidade Barras PV
num_n_PQ    = length(Index_n_PQ);                               % Quantidade Barras PQ
num_n_Sl    = num_barras - (num_n_PV + num_n_PQ);               % Quantidade Barras Sl



%||  Criar barras virtuais em nodos PV ou slack que seja precisado  ||
num_n_novos = 0;
for i = 1 : num_n_PV
    if( Dados_barras( Index_n_PV(i), 6 ) ~= 0)                  % Nodos PV com carga Q
        novo_n                        = aux_novo_n;
                                                         %num   tip Pg  Qg  Pc   Qc
        Dados_barras                  = [Dados_barras ; novo_n , 3 , 0 , 0 , 0 , Dados_barras( Index_n_PV(i), 6 ), 1 ,0 ];
        Dados_barras(Index_n_PV(i),6) = 0;
        aux_novo_n                    = aux_novo_n + 1;
        num_n_novos                   = num_n_novos + 1;
        num_n_PQ                      = num_n_PQ + 1;
        num_barras                    = num_barras + 1;
        num_linhas                    = num_linhas + 1;
                                                         %ni                            nj      R     X        Y/2 
        Dados_linhas                  = [Dados_linhas; Dados_barras(Index_n_PV(i),1) , novo_n , 0 , 0.00001 , 0];
    end
end

for i = 1 : num_n_Sl
    if( (Dados_barras( Index_n_Sl(i), 5 ) ~= 0) || (Dados_barras( Index_n_Sl(i), 6 ) ~= 0) )    % Nodos Sl com carga P ou Q
        novo_n                        = aux_novo_n;
                                                         %num   tip Pg  Qg  Pc   Qc
        Dados_barras                  = [Dados_barras ; novo_n , 3 , 0 , 0 , Dados_barras( Index_n_Sl(i), 5 ) , Dados_barras( Index_n_Sl(i), 6 ), 1 ,0 ];
        Dados_barras(Index_n_Sl(i),5) = 0;
        Dados_barras(Index_n_Sl(i),6) = 0;
        aux_novo_n                    = aux_novo_n + 1;
        num_n_novos                   = num_n_novos + 1;
        num_n_PQ                      = num_n_PQ + 1;
        num_barras                    = num_barras + 1;
        num_linhas                    = num_linhas + 1;
                                                         %ni                            nj      R     X        Y/2 
        Dados_linhas                  = [Dados_linhas; Dados_barras(Index_n_Sl(i),1) , novo_n , 0 , 0.00001 , 0];
    end
end
clear aux_novo_n;
clear novo_n;
%---------------------------------------------------------------------------------------------------------------------

%||  Renumeração processo aux||
%==============================

pos2        = num_n_PQ + num_n_PV;
pos3        = num_n_PQ;

vetor_aux   = unique([Dados_linhas(:,1) ; Dados_linhas(:,2)].');     % Coloca em ordem ascendente e sem repetir; Transpuesto
vetor_aux2  = zeros(1, num_barras);                                 % Vetor formado para realizar o ordenamento colocando primeiro as barras PQ, depois PV e de último slack

for i = 1 : num_barras
    for j = 1 : num_barras
        
          if( Dados_barras(i,1) == vetor_aux(j) )               %Posição "i" do vetor aux2 tem o número e posição da barra da posição "i" do vetor aux
              
              if( Dados_barras(i,2) == 1)
                  vetor_aux2(j) = num_barras;
                  break
              elseif( Dados_barras(i,2) == 2 )
                  vetor_aux2(j) = pos2;
                  pos2           = pos2-1;
                  break
              elseif( Dados_barras(i,2) == 3 )
                  vetor_aux2(j) = pos3;
                  pos3           = pos3-1;
                  break
              end
         
          end
    end       
end
clear pos2;
clear pos3


%||  Renumeração ||
%==================

for i = 1 : num_linhas
    flag1 = 1;
    flag2 = 1;
    flag3 = 1;
    for j = 1 : num_barras
        
        if( flag1 && Dados_linhas(i,1) == vetor_aux(j))           
             Dados_linhas(i,1) = vetor_aux2(j);
             flag1 = 0;
             if (flag1==0 && flag2==0 && flag3==0)
                 break;
             end
        end
        
        if( flag2 && Dados_linhas(i,2) == vetor_aux(j))           
             Dados_linhas(i,2) = vetor_aux2(j);
             flag2 = 0;
             if (flag1==0 && flag2==0 && flag3==0)
                 break;
             end
        end
        
        if( flag3 && (i <= num_barras) && (Dados_barras(i,1) == vetor_aux(j)) )  
             Dados_barras(i,1) = vetor_aux2(j);
             flag3 = 0;
             if (flag1==0 && flag2==0 && flag3==0)
                 break;
             end
        end
    end
end
clear vetor_aux;
clear flag1;
clear flag2;
clear flag3;

end

