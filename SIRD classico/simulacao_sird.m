% === SIMULAÇÃO DO MODELO SIRD - BAURU (400 dias) ===

% Parâmetros exatos fornecidos
beta   = 0.000000294858301579693;    % Taxa de contato
lambda = 0.0345122943418598;         % Taxa de recuperação
gamma  = 0.00835692416616626;        % Taxa de mortalidade


% População total
N = 379300;

% Condições iniciais
S0 = 379298;
I0 = 2;
R0 = 0;
D0 = 0;
Y0 = [S0; I0; R0; D0];

% Intervalo de tempo da simulação (400 dias)
tspan = [0 400];

% Resolver o sistema com ode45
[t, Y] = ode45(@(t, Y) sird_model(t, Y, beta, lambda, gamma), tspan, Y0);

% Separar as soluções
S = Y(:,1);
I = Y(:,2);
R = Y(:,3);
D = Y(:,4);

% Plot completo com todas as curvas
figure;
plot(t, S, 'b', 'LineWidth', 2); hold on;
plot(t, I, 'r', 'LineWidth', 2);
plot(t, R, 'g', 'LineWidth', 2);
plot(t, D, 'k', 'LineWidth', 2);
xlabel('Dias desde 30/03/2020');
ylabel('Número de Pessoas');
legend('Suscetíveis', 'Infectados', 'Recuperados', 'Óbitos', 'Location', 'northeast');
title('Simulação do Modelo SIRD - Bauru (400 dias)');
grid on;
