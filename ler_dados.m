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
%Função para Ler número de linhas e barras e Pbase, Vbase.

function [P_base, V_base, num_linhas, num_barras, Dados_linhas, Dados_barras] = ler_dados(nome_sis)
%||  Ler dados ||
%================
Valrs_base = xlsread(nome_sis, 'Hoja1', 'B1:B3');                   %Por enquanto não ingresar ao programa barras slack com cergas 
P_base     = Valrs_base(1);
V_base     = Valrs_base(2);
num_linhas = uint8(Valrs_base(3));                                  %Convertir a valor inteiro sem signo
flag_pu    = xlsread(nome_sis, 'Hoja1', 'E1:E1');
%num_barras = uint8(Valrs_base(4));                                 %Convertir a valor inteiro sem signo
clear Valrs_base;

Dados_linhas = xlsread(nome_sis, 'Hoja1', strcat( 'A8:E',num2str(num_linhas+7) ));

vetor_aux  = unique([Dados_linhas(:,1) ; Dados_linhas(:,2)].');     % Coloca em ordem ascendente e sem repetir; Transpuesto
num_barras = length(vetor_aux);                                     % Conhecer num de barras

%---------------------------------------------------------------------------------------------------------------------
Dados_barras = xlsread(nome_sis, 'Hoja1', strcat( 'A',num2str(num_linhas+10),':H',num2str(num_linhas+10+num_barras) ));
%||  Colocar Potências em p.u. ||
if(~flag_pu)
    Dados_barras(:, [3,4,5,6]) = (1/P_base)*Dados_barras(:, [3,4,5,6]);
end
clear flag_pu;
%---------------------------------------------------------------------------------------------------------------------




end
