%  <CÁLCULO DO FLUXO DE POTÊNCIA - DC LOAD FLOW V1.0. 
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

%Cálculo_Fluxo_de_Potência_linearizado
%   Sebastián de Jesús Manrique Machado
%   Estudante_Mestrado Em Engenharia Elétrica
%   UEL - 2014.

nome_sis   = input('Qual é o nome do sistema para executar FP?: ','s'); 
tol        = 0.0005;
cont       = 1;

disp('Inicio Cálculo Fluxo Linearizado DC');
disp(datestr(now));

%||  Ler dados ||
%================
[P_base, V_base, num_linhas_or, num_barras_or, Dados_linhas_or, Dados_barras_or] = ler_dados(nome_sis);

[ Dados_linhas, Dados_barras, index_original, num_barras, num_linhas, num_n_PV, num_n_PQ, num_n_novos ] = renumerar( Dados_linhas_or, Dados_barras_or, num_linhas_or, num_barras_or );

nodo_i      = Dados_linhas(:,1);
nodo_j      = Dados_linhas(:,2);

%|| Ordenamento de dados ||                                 Num da barra = posição no vetor
%==========================
Barra       = zeros(num_barras, 1);
Tipo_barra  = zeros(num_barras, 1);                        %1=Slack; 2=PV; 3=PQ
P_dado      = zeros(num_barras, 1);                        %Pgen-Pcarga 
P_dado_p    = zeros(num_barras, 1);                        %Pgen-Pcarga-Pperdas
v_barras    = ones(num_barras, 1);                          % No fluxo DC todas as tensôes são consideradas 1 p.u.
delt_barras = zeros(num_barras, 1);

for i = 1 : num_barras
    for j = 1 : num_barras
        if( Dados_barras(j,1) == i )
            Barra(i)       = i;                                        
            Tipo_barra(i)  = Dados_barras(j,2);                        %1=Slack; 2=PV; 3=PQ
            P_dado(i)      = Dados_barras(j,3) - Dados_barras(j,5);    %Pgen-Pcarga
            P_dado_p(i)    = P_dado(i);
            %v_barras(i)    = Dados_barras(j,7);
            delt_barras(i) = Dados_barras(j,8);
        end
    end
end


%||  Cálculo Y Barras ||
%=======================
disp('Cálculo Ybarras')

%Levando em conta a resistência
Zs_linha      = Dados_linhas(:,3) + 1i*Dados_linhas(:,4);
Ys_linha      = 1./Zs_linha;
b_shunt       = 1i*Dados_linhas(:,5);
[Y_barras]    = calculo_Yb(Ys_linha, b_shunt, num_barras, num_linhas, nodo_i, nodo_j);
B_barras_p    = imag( Y_barras(1:num_n_PV+num_n_PQ, 1:num_n_PV+num_n_PQ) );



delt_barras(1:num_n_PV+num_n_PQ , 1:1) = -B_barras_p\P_dado(1:num_n_PV+num_n_PQ , 1:1);

flag_p      = input('Deseja fazer fluxo estimando perdas? (1=sim; cualquer outro número = não): ');          

if(flag_p)
    cont          = 2;
    G_barras_p    = real(Y_barras);
    Perdas_ij     = zeros(num_barras, num_barras);
    for i = 1 : num_linhas
    Perdas_ij(nodo_i(i), nodo_j(i)) = abs(G_barras_p(nodo_i(i), nodo_j(i))) * (delt_barras( Barra==nodo_i(i) )-delt_barras( Barra==nodo_j(i) ))^2;
    Perdas_ij(nodo_j(i), nodo_i(i)) = Perdas_ij(nodo_i(i), nodo_j(i));
    P_dado_p(nodo_i(i))             = P_dado_p(nodo_i(i)) - ( Perdas_ij(nodo_i(i), nodo_j(i)) / 2 );
    P_dado_p(nodo_j(i))             = P_dado_p(nodo_j(i)) - ( Perdas_ij(nodo_j(i), nodo_i(i)) / 2 );
    end
    delt_barras(1:num_n_PV+num_n_PQ , 1:1) = -B_barras_p\P_dado(1:num_n_PV+num_n_PQ , 1:1);
end


V_complejo  = v_barras.* cos(delt_barras) + 1i * v_barras.* sin(delt_barras);
I_inj       = Y_barras * V_complejo;



%||  Cálculo subsistema 2 S_Linhas  ||
%=====================================
S_ij      = zeros(num_barras, num_barras);
I_ij      = zeros(num_barras, num_barras);
%||  Cálculo I e S pelas linhas ||
for i = 1 : num_linhas
    I_ij(nodo_i(i), nodo_j(i)) = (- V_complejo(nodo_i(i)) + V_complejo(nodo_j(i)) ) * Y_barras(nodo_i(i), nodo_j(i)) + V_complejo(nodo_i(i)) * b_shunt(i);
    I_ij(nodo_j(i), nodo_i(i)) = (- V_complejo(nodo_j(i)) + V_complejo(nodo_i(i)) ) * Y_barras(nodo_j(i), nodo_i(i)) + V_complejo(nodo_j(i)) * b_shunt(i);
    
    S_ij(nodo_i(i), nodo_j(i)) = (delt_barras( Barra==nodo_i(i) )-delt_barras( Barra==nodo_j(i) )) / abs(imag(Zs_linha(i)));
    S_ij(nodo_j(i), nodo_i(i)) = -S_ij(nodo_i(i), nodo_j(i));
    if(S_ij(nodo_i(i), nodo_j(i)) > 0)
        S_ij(nodo_j(i), nodo_i(i)) = S_ij(nodo_j(i), nodo_i(i)) + Perdas_ij(nodo_j(i), nodo_i(i));
    else
        S_ij(nodo_i(i), nodo_j(i)) = S_ij(nodo_i(i), nodo_j(i)) + Perdas_ij(nodo_i(i), nodo_j(i));
    end
    if((nodo_i(i) == num_barras) || (nodo_j(i) == num_barras))
        P_dado(num_barras) = P_dado(num_barras) + S_ij(nodo_i(i), nodo_j(i));
    end
end

for i = 1 : num_barras
    I_ij(i , i) = I_inj(i);
    S_ij(i , i) = P_dado(i);
end





%||  Imprimir arquivo de texto  ||
%=================================
imprimir_res( index_original, num_barras, num_n_novos, num_linhas, nodo_i, nodo_j, P_dado, P_dado-P_dado, v_barras, delt_barras, S_ij, Perdas_ij, I_ij, P_base, V_base, cont, nome_sis, '4linear',Dados_linhas_or, Dados_barras_or, num_n_PQ )


disp('Terminou Cálculo Fluxo');
disp(datestr(now));
