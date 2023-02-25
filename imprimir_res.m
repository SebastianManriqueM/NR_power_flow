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
%Função para imprimir arquivo de texto.

function imprimir_res( index_original, num_barras, num_n_novos , num_linhas, nodo_i, nodo_j, P_calc, Q_calc, v_barras, delt_barras, S_ij, Loses_ij, I_ij, P_base, V_base, cont, nome_sis, metodo, Dados_linhas_or, Dados_barras_or, num_n_PQ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%------------------------------------------------
%|||||||    Impresão dados de Entrada    |||||||
%------------------------------------------------
tex=fopen(['Resultado_fluxo_', nome_sis,'-', metodo, '.txt'],'w');                  %Abro arquivo de texto
fprintf(tex, '========================================\r\n');
fprintf(tex, '|||||   DADOS BÁSICOS DE ENTRADA   |||||\r\n');
fprintf(tex, '========================================\r\n\r\n');
fprintf(tex,strcat('\r\tData: ', datestr(now), '\r\n\r\n'));
fprintf(tex,strcat('\r\tNome do sistema: ', nome_sis,'\r\n') );
fprintf(tex,strcat('\r\tMétodo Usado: ', metodo,'\r\n') );
fprintf(tex,strcat('\r\tPotência base= ', num2str(P_base), 'MW', '\r\t', 'Tensão Base= ', num2str(V_base),'kV', '\r\n') );
fprintf(tex,strcat('\r\tBarras= ', num2str(num_barras-num_n_novos), '\r\t', '\r\n') );
fprintf(tex,strcat('\r\tNumero de Linhas: ', num2str(num_linhas-num_n_novos),'\r\n\r\n') );
fprintf(tex,strcat( '\r\tITERAÇÕES =  ', num2str(cont), '\r\n\r\n\r\n' ));
%:::::::::::::::::::::::::::::::::::::::::::::::::::

%------------------------------------------------
%|||||||    Impresão Resultados Linhas    |||||||
%------------------------------------------------
S_ij     = P_base * S_ij;
Loses_ij = P_base * Loses_ij;

fprintf(tex, '=================================\r\n');
fprintf(tex, '|||||   RESULTADOS LINHAS   |||||\r\n');
fprintf(tex, '=================================\r\n\r\n');
fprintf(tex, 'Bus i\r\tBus j\r\t\r\tPij [MW]\r\tQij [Mvars]\r\tPji [MW] \r\tQji [Mvars]\r\t|Iij| [pu]\r\t<Iij [º]\r\tPerdas [MW]\r\tPerdas [MVars]\r\n\r\n' );

for i = 1 : num_linhas-num_n_novos
    fprintf(tex,'%3.1f\r\t %3.1f\r\t %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\n',Dados_linhas_or(i,1), Dados_linhas_or(i,2), real( S_ij(nodo_i(i),nodo_j(i)) ), imag( S_ij(nodo_i(i),nodo_j(i)) ), real( S_ij(nodo_j(i),nodo_i(i)) ), imag( S_ij(nodo_j(i),nodo_i(i)) ), abs(I_ij(nodo_i(i),nodo_j(i))), radtodeg(angle(I_ij(nodo_i(i),nodo_j(i)))), real( Loses_ij(nodo_i(i),nodo_j(i)) ), imag( Loses_ij(nodo_i(i),nodo_j(i)) ) );
end

%------------------------------------------------
%|||||||    Impresão Resultados Barras    |||||||
%------------------------------------------------
P_calc = P_base * P_calc;
Q_calc = P_base * Q_calc;
Dados_barras_or(:, [3,4,5,6]) = (P_base)*Dados_barras_or(:, [3,4,5,6]);

fprintf(tex, '\r\n\r\n\r\n=================================\r\n');
fprintf(tex, '|||||   RESULTADOS BARRAS   |||||\r\n');
fprintf(tex, '=================================\r\n\r\n');
fprintf(tex, 'Bus i\r\t\r\tV [p.u.]\r\tdelta [º]\r\tPg [MW]  \r\tQg [Mvars]   \r\tPl [MW]  \r\tQl [Mvars]   \r\tPc [MW]  \r\tQc [Mvars]\r\n\r\n' );
cont = 2;   % nodos slack + 1!!!.
for i = num_n_novos+1 : num_barras
     a = find(index_original==i);
     if (i == num_barras)                   %Dados Barra Slack
        fprintf(tex,'%3.1f\r\t %12.3f\r\t  %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\n', Dados_barras_or(a,1), v_barras(i), radtodeg(delt_barras(i)), P_calc(i), Q_calc(i), Dados_barras_or(a,5), Dados_barras_or(a,6), P_calc(i), Q_calc(i) );
     elseif(i > num_n_PQ)                   %Dados Barras PV
        fprintf(tex,'%3.1f\r\t %12.3f\r\t  %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\n', Dados_barras_or(a,1), v_barras(i), radtodeg(delt_barras(i)), Dados_barras_or(a,3), Q_calc(i), Dados_barras_or(a,5), Dados_barras_or(a,6), P_calc(i), Q_calc(i) );         
        cont = cont + 1;
     else
        fprintf(tex,'%3.1f\r\t %12.3f\r\t  %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\t   %12.3f\r\n', Dados_barras_or(a,1), v_barras(i), radtodeg(delt_barras(i)),Dados_barras_or(a,3), Dados_barras_or(a,4), Dados_barras_or(a,5), Dados_barras_or(a,6), P_calc(i), Q_calc(i) );
        
      end
end

fprintf(tex, '\r\n\r\n**As barras do sistema forãm renumeradas de tal forma que a última barra é a barra slack, e as primeiras barras são as PQ.\r\nCopyright (C) <2014>  <Sebastián de Jesús Manrique Machado>   <e-mail:sebajmanrique747@gmail.com>');