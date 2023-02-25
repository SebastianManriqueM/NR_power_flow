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
%Função para achar Jacobiano.

function [ J_11, J_12, J_21, J_22 ] = jacobiano( v_barras, delt_barras, G_barras, B_barras, P_calc, Q_calc, num_n_PV, num_n_PQ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

J_11 = zeros(num_n_PV + num_n_PQ , num_n_PV + num_n_PQ);

%||  Derivada de P respeito a Theta ||      H
%=====================================
% Em barras PV e PQ tem-se P dado e a variável Theta é desconhecida sempre
for i = 1 : num_n_PV + num_n_PQ
    
    for j = 1 : num_n_PV + num_n_PQ
        
        if (i == j)
            J_11(i , j) = ( -(v_barras(i))^2 ) * B_barras(i , j) - Q_calc(i);
        elseif(i > j)
            J_11(i , j) =  v_barras(i) * v_barras(j) * ( (G_barras(i,j)*sin(delt_barras(i)-delt_barras(j)) ) - (B_barras(i,j)*cos(delt_barras(i)-delt_barras(j)) )  );
            J_11(j , i) =  v_barras(i) * v_barras(j) * ( (G_barras(j,i)*sin(delt_barras(j)-delt_barras(i)) ) - (B_barras(j,i)*cos(delt_barras(j)-delt_barras(i)) )  );
        end
        
    end
    
end


%||  Derivada de P respeito a V ||      N
%=================================
% Em barras PV e PQ tem-se P dado e a variável V é desconhecida nas PQ
J_12 = zeros(num_n_PV + num_n_PQ , num_n_PQ);

for i = 1 : num_n_PV + num_n_PQ
    
    for j = 1 : num_n_PQ
        
        if (i == j)
            J_12(i , j) = ( v_barras(i) * G_barras(i , j)) + (1/v_barras(i))*P_calc(i);
        else 
            J_12(i , j) =  v_barras(i) * ( (G_barras(i,j)*cos(delt_barras(i)-delt_barras(j)) ) + (B_barras(i,j)*sin(delt_barras(i)-delt_barras(j)) )  );
        end
        
    end
    
end


%||  Derivada de Q respeito a Theta ||
%=====================================
J_21 = zeros(num_n_PQ , num_n_PV + num_n_PQ);

for i = 1 : num_n_PQ
    
    for j = 1 : num_n_PV + num_n_PQ
        
        if (i == j)
            J_21(i , j) = P_calc(i)-(v_barras(i)^2*G_barras(i,j));
        else  
            J_21(i , j) =  -v_barras(i) * v_barras(j) *( (G_barras(i,j)*cos(delt_barras(i)-delt_barras(j)) ) + (B_barras(i,j)*sin(delt_barras(i)-delt_barras(j)) )  );
        end
        
    end
    
end



%||  Derivada de Q respeito a V ||
%=================================
J_22 = zeros(num_n_PQ , num_n_PQ);

for i = 1 : num_n_PQ
    
    for j = 1 : num_n_PQ
        
        if (i == j)
            J_22(i , j) = ( -(v_barras(i)) * B_barras(i , j))  + (1/v_barras(i))*Q_calc(i);
        elseif(i > j)
            J_22(i , j) =  v_barras(i) * ( (G_barras(i,j)*sin(delt_barras(i)-delt_barras(j)) ) - (B_barras(i,j)*cos(delt_barras(i)-delt_barras(j)) )  );
            J_22(j , i) =  v_barras(j) * ( (G_barras(j,i)*sin(delt_barras(j)-delt_barras(i)) ) - (B_barras(j,i)*cos(delt_barras(j)-delt_barras(i)) )  );
        end
        
    end
    
end


end

