%% ORDENS DAS DERIVADAS
alpha1 = 0.99;   % Ordem derivada S
alpha2 = 0.99;   % Ordem derivada I
alpha3 = 0.99;   % Ordem derivada R
alpha4 = 0.99;   % Ordem derivada D

alpha = [alpha1, alpha2, alpha3, alpha4 ]; % Vetor ordem das derivadas

%% PARÂMETROS
beta = 0.000000294858301579693 ; % Parâmetro beta -   Taxa de Contato
lambda = 0.0345122943418598;     % Parâmetro lambda - Taxa de Recuperação
gama = 0.00835692416616626 ;     % Parâmetro gamma -  Taxa de Mortalidade

param = [ beta , lambda  , gama  ] ;  % Vetor de parâmetros

h=0.02;      % Passo
t0=0;        % Tempo inicial           
tf=400;      % Tempo final                   
t=t0:h:tf;   % Intervalo de tempo

%% CONDIÇÕES INICIAS
S0= 379298;     % Numero Inicial de Suscetiveis
I0 = 2;         % Numero Inicial de Infectados
R0 = 0;         % Numero Inicial de Recuperados
D0 = 0;         % Numero Inicial de Mortos
y0 = [S0 ; I0 ; R0 ; D0]; % Vetor condição inicial
dS0 = -5.89716603159386e-07;
dI0 = -0.08573784729944896;
dR0 =  0.0690245886837196;
dD0 =  0.01671384833233252;
dy0 = [dS0; dI0; dR0; dD0];
Yi = [y0, dy0];  % Onde y0 = [S0; I0; R0; D0] e dy0 =  [dS0; dI0; dR0; dD0]


%% FUNÇÕES
f_fun = @(t,y,param) [ -(param(1)^alpha1)*y(1)*y(2); % Equação dos Suscestíveis
                     (param(1)^alpha2)*y(1)*y(2) - (param(2)^alpha2)*y(2)- (param(3)^alpha2)*y(2); % Equação dos Infectados
                     (param(2)^alpha3)*y(2); % Equação dos Recuperados
                     (param(3)^alpha4)*y(2) ]; % Equação dos Mortos
J_fun = @(t, y, param) [  
    -(param(1)^alpha1)*y(2),                  -(param(1)^alpha1)*y(1),                          0,                          0; 
     (param(1)^alpha2)*y(2),   (param(1)^alpha2)*y(1) - (param(2)^alpha2) - (param(3)^alpha2),  0,                          0; 
     0,                                      (param(2)^alpha3),                                 0,                          0; 
     0,                                      (param(3)^alpha3),                                 0,                          0  ]; 

%% ROTINAS
%[f, y] = fde_pi1_ex(alpha,f_fun,t0,tf,Yi,h,param);       % explicit product-integration rule of rectangular type with convergence order equal to 1.      
[f, y] = fde_pi12_pc(alpha,f_fun,t0,tf,y0,h,param);      % is solved by means of a couple of product-integration rules of rectangular and trapezoidal type implemented in a predictor-corrector framework.  
%[f, y] = fde_pi1_im(alpha,f_fun,J_fun,t0,tf,Yi,h,param);  % implicit product-integration rule of trapezoidal type with convergence order equal to 1 (under suitable assumptions of regularity of the exact solution).
%[f, y] = fde_pi2_im(alpha,f_fun,J_fun,t0,tf,Yi,h,param); % implicit product-integration rule of trapezoidal type with convergence order equal to 2 (under suitable assumptions of regularity of the exact solution).
%[tc, y] = fde_pi1_ex(alpha2,f,t02,tf2,N02,h2,w2);
% SAÍDAS
S= y(1,:);     % Resultados Suscetíveis
I= y(2,:);  % Resultados Infectados
R= y(3,:);  % Resultados Recuperados
D= y(4,:);  % Resultados Mortos

%% GRÁFICOS
figure(1)
plot(t,S,'b',t,I,'m',t,R,'g',t,D,'r','LineWidth', 1.5);
grid on;
%ylim([0 15000]);
title('Modelo SIRD Fracionário com \alpha = 1.02');
xlabel('Tempo (em dias)');                                           
ylabel('População (nº de indivíduos)');
hold on;
legend('Suscetíveis','Infectados','Recuperados', 'Mortos');
hold on;

% plot(t, I, 'Color', colors(k,:), 'LineWidth', 1.5, 'DisplayName', ['\alpha = ', num2str(alpha_values(k))]);

%%



















