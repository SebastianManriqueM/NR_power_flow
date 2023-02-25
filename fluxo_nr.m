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

nome_sis   = input('Qual é o nome do sistema para executar FP?: ','s'); %Datos_ej1.xlsx 
tol        = 0.0005;
disp('Inicio Cálculo Fluxo Newton Raphson');
disp(datestr(now));

%||  Ler dados ||
%================
[P_base, V_base, num_linhas_or, num_barras_or, Dados_linhas_or, Dados_barras_or] = ler_dados(nome_sis);

[ Dados_linhas, Dados_barras, index_original, num_barras, num_linhas, num_n_PV, num_n_PQ, num_n_novos ] = renumerar( Dados_linhas_or, Dados_barras_or, num_linhas_or, num_barras_or );

nodo_i      = Dados_linhas(:,1);
nodo_j      = Dados_linhas(:,2);
Zs_linha    = Dados_linhas(:,3) + 1i*Dados_linhas(:,4);
Ys_linha    = 1./Zs_linha;                              %Admitância Serie da Linha
b_shunt     = 1i*Dados_linhas(:,5);

tol_P = zeros(num_n_PV + num_n_PQ, 1);
for i = 1 : num_n_PV + num_n_PQ
    tol_P(i) = tol;
end

tol_Q = zeros(num_n_PQ, 1);
for i = 1 : num_n_PQ
    tol_Q(i) = tol;
end

%|| Ordenamento de dados ||                                 Num da barra = posição no vetor
%==========================
Barra       = zeros(num_barras, 1);
Tipo_barra  = zeros(num_barras, 1);                        %1=Slack; 2=PV; 3=PQ
P_dado      = zeros(num_barras, 1);                        %Pgen-Pcarga
Q_dado      = zeros(num_barras, 1);                        %Qgen-Qcarga
v_barras    = zeros(num_barras, 1);
delt_barras = zeros(num_barras, 1);

for i = 1 : num_barras
    for j = 1 : num_barras
        if( Dados_barras(j,1) == i )
            Barra(i)       = i;                                        
            Tipo_barra(i)  = Dados_barras(j,2);                        %1=Slack; 2=PV; 3=PQ
            P_dado(i)      = Dados_barras(j,3) - Dados_barras(j,5);    %Pcarga-Pgen
            Q_dado(i)      = Dados_barras(j,4) - Dados_barras(j,6);    %Qgen-Qcarga
            v_barras(i)    = Dados_barras(j,7);
            delt_barras(i) = Dados_barras(j,8);
            break;
        end
    end
end


%||  Cálculo Y Barras ||
%=======================
disp('Cálculo Ybarras')
[Y_barras] = calculo_Yb(Ys_linha, b_shunt, num_barras, num_linhas, nodo_i, nodo_j);
G_barras   = real(Y_barras);
B_barras   = imag(Y_barras);

Delta_theta = zeros(num_n_PV + num_n_PQ, 1);
Delta_V     = zeros(num_n_PQ, 1);
Delta_x     = [Delta_theta; Delta_V];


%||  ITERAÇÕES ||
%================
flag = 1;
cont = 1;
while(flag)
    disp( strcat('Iteração ',num2str(cont)) )
    disp('----------')
    %||  1  Cálculo do valor da função no ponto atual ||
    V_complejo  = v_barras.* cos(delt_barras) + 1i * v_barras.* sin(delt_barras);   % Não é produto matrizial. Produto termo a termo
    I_inj  = Y_barras * V_complejo;
    S_calc = V_complejo.* conj(I_inj);                                              % Não é produto matrizial. Produto termo a termo
    P_calc = real(S_calc);
    Q_calc = imag(S_calc);
    
    Delta_P = P_dado - P_calc;
    Delta_Q = Q_dado - Q_calc;  
    
    %||  2  Cumpre com a Tolerância? ||
    if( all(abs(Delta_P(1:num_n_PV+num_n_PQ , 1)) <= tol_P) && all(abs(Delta_Q(1:num_n_PQ , 1)) <= tol_Q) )
        flag = 0;
        break
    end
    
    
    %||  3  Cálculo Jacobiano ||
    [ J_11, J_12, J_21, J_22 ]   = jacobiano( v_barras, delt_barras, G_barras, B_barras, P_calc, Q_calc, num_n_PV, num_n_PQ );
    J          = [J_11 J_12; J_21 J_22];
    
    %||  4  Cálculo Delta theta e Delta V ||
   Delta_x     = J\[Delta_P(1:num_n_PV+num_n_PQ , 1); Delta_Q(1:num_n_PQ , 1)];     %inv de J = J\
   Delta_theta = Delta_x(1:num_n_PV+num_n_PQ,1);
   Delta_V     = Delta_x(num_barras:num_n_PV+2*num_n_PQ,1);
   
   
    %||  5  Novo vetor de variáveis ||
   v_barras(1:num_n_PQ , 1)             = v_barras(1:num_n_PQ,1) + Delta_V;
   delt_barras(1:num_n_PV+num_n_PQ , 1) = delt_barras(1:num_n_PV+num_n_PQ , 1) + Delta_theta;
   
   cont = cont + 1;
end 


%||  Cálculo subsistema 2 S_Linhas  ||
%=====================================
[S_ij, I_ij, Perdas_ij] = calculo_linhas(num_barras, num_linhas, nodo_i, nodo_j, Y_barras, b_shunt, I_inj, V_complejo, S_calc);


%||  Imprimir arquivo de texto  ||
%=================================
imprimir_res( index_original, num_barras, num_n_novos, num_linhas, nodo_i, nodo_j, P_calc, Q_calc, v_barras, delt_barras, S_ij, Perdas_ij, I_ij, P_base, V_base, cont, nome_sis, '1nr', Dados_linhas_or, Dados_barras_or, num_n_PQ )


disp('Terminou Cálculo Fluxo');
disp(datestr(now));
