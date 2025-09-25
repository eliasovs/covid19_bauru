% === COMPARAÇÃO ENTRE ESTRATÉGIAS 2, 3 E 4 ===

% --- Dados reais (1:62)
t = 0:61;  % Dias
casos_reais = casos_acumulados(1:62);
obitos_reais = obitos_acumulados(1:62);

% === Estratégia 2 ===
beta2   = 0.000001560969841;
lambda2 = 0.174349279428620;
gamma2  = 0.014894814476100;
alpha2  = 1.10835;

f2 = @(t, y, p) [
    -(p(1)^alpha2) * y(1) * y(2);
     (p(1)^alpha2) * y(1) * y(2) - (p(2)^alpha2) * y(2) - (p(3)^alpha2) * y(2);
     (p(2)^alpha2) * y(2);
     (p(3)^alpha2) * y(2)
];

y0 = [379298; 2; 0; 0];
dy0_2 = [
    -beta2 * y0(1) * y0(2);
     beta2 * y0(1) * y0(2) - lambda2 * y0(2) - gamma2 * y0(2);
     lambda2 * y0(2);
     gamma2 * y0(2)
];
y_ini2 = [y0, dy0_2];
[t_sim, y2] = fde_pi12_pc([alpha2, alpha2, alpha2, alpha2], f2, 0, 61, y_ini2, 1, [beta2, lambda2, gamma2]);
casos_modelo_estrategia2  = y2(2,:) + y2(3,:) + y2(4,:);
obitos_modelo_estrategia2 = y2(4,:);

% === Estratégia 3 ===
beta3   = 0.000000609831270;
lambda3 = 0.098628811273860;
gamma3  = 0.012132630254103;
alpha13 = 0.35225;
alpha23 = 1.00493;
alpha33 = 1.07332;
alpha43 = 1.00000;

f3 = @(t, y, p) [
    -(p(1)^alpha13) * y(1) * y(2);
     (p(1)^alpha23) * y(1) * y(2) - (p(2)^alpha23) * y(2) - (p(3)^alpha23) * y(2);
     (p(2)^alpha43) * y(2);
     (p(3)^alpha33) * y(2)
];

dy0_3 = [
    -beta3 * y0(1) * y0(2);
     beta3 * y0(1) * y0(2) - lambda3 * y0(2) - gamma3 * y0(2);
     lambda3 * y0(2);
     gamma3 * y0(2)
];
y_ini3 = [y0, dy0_3];
[~, y3] = fde_pi12_pc([alpha13, alpha23, alpha33, alpha43], f3, 0, 61, y_ini3, 1, [beta3, lambda3, gamma3]);
casos_modelo_estrategia3  = y3(2,:) + y3(3,:) + y3(4,:);
obitos_modelo_estrategia3 = y3(4,:);

% === Estratégia 4 ===
beta4   = 0.000001614274829;
lambda4 = 0.095357881104545;
gamma4  = 0.049999999922082;
alpha14 = 1.00100;
alpha24 = 1.20241;
alpha34 = 1.50000;
alpha44 = 1.05527;

f4 = @(t, y, p) [
    -(p(1)^alpha14) * y(1) * y(2);
     (p(1)^alpha24) * y(1) * y(2) - (p(2)^alpha24) * y(2) - (p(3)^alpha24) * y(2);
     (p(2)^alpha44) * y(2);
     (p(3)^alpha34) * y(2)
];

dy0_4 = [
    -beta4 * y0(1) * y0(2);
     beta4 * y0(1) * y0(2) - lambda4 * y0(2) - gamma4 * y0(2);
     lambda4 * y0(2);
     gamma4 * y0(2)
];
y_ini4 = [y0, dy0_4];
[~, y4] = fde_pi12_pc([alpha14, alpha24, alpha34, alpha44], f4, 0, 61, y_ini4, 1, [beta4, lambda4, gamma4]);
casos_modelo_estrategia4  = y4(2,:) + y4(3,:) + y4(4,:);
obitos_modelo_estrategia4 = y4(4,:);

% === GRÁFICO: CASOS ACUMULADOS ===
figure;
plot(t, casos_reais, 'ko', 'DisplayName', 'Dados reais'); 
hold on;
plot(t, casos_modelo_estrategia2, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2');
plot(t, casos_modelo_estrategia3, 'c-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 3');
plot(t, casos_modelo_estrategia4, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 4');
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11); ylabel('Casos acumulados', 'FontSize', 11);
title('Casos acumulados de Covid-19 em Bauru: Comparação entre Estratégias 2, 3 e 4', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12); 
grid on;



%% === GRÁFICO: CASOS ACUMULADOS === COMPARANDO COM CLASSICO
figure;
hold on; % <-- aqui!
%plot(t, casos_classico, 'k--', 'LineWidth', 2, 'DisplayName', 'Curva estimada de casos acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
%plot(t, casos_101, 'c-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
%plot(t, casos_modelo, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
%plot(t, casos_modelo_estrategia2, 'b-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2');
plot(t, casos_modelo_estrategia4, 'g-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 4');
plot(t, casos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de casos acumulados de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11);
ylabel('Casos acumulados', 'FontSize', 11);
title('Casos acumulados de Covid-19 em Bauru', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12); 
ylim([0 400]);  % <-- Aqui você fixa os limites do eixo Y
grid on;

%% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
hold on;
plot(t, obitos_modelo_estrategia2, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2');
%plot(t, obitos_modelo_estrategia3, 'y-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 3');
plot(t, obitos_modelo_estrategia4, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 4');
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11); ylabel('Casos acumulados', 'FontSize', 11);
title('Mortes Acumuladas de Covid-19 em Bauru: Comparação entre Estratégias 2, 3 e 4', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12); 
grid on;


%% === GRÁFICO: MORTES ACUMULADAS === COMPARANDO COM CLASSICO

figure;
hold on;
%plot(t, obitos_classico, 'k--', 'LineWidth', 2, 'DisplayName', 'Curva estimada de mortes acumulados de Covid-19 em Bauru: modelo clássico (\alpha = 1)');
%plot(t, obitos_modelo_estrategia2, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2');
%plot(t, obitos_modelo_estrategia3, 'y-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 3');
plot(t, obitos_modelo_estrategia4, 'm-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 4');
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumuladas de Covid-19 em Bauru');
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11); ylabel('Casos acumulados', 'FontSize', 11);
title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12); 
ylim([0 25]);  % <-- Aqui você fixa os limites do eixo Y
grid on;



%%
% === GRÁFICO: ÓBITOS ACUMULADOS ===
figure;
hold on;

%plot(t, obitos_101, 'y-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 1 \alpha = 1.01');
plot(t, obitos_modelo, 'r-', 'LineWidth', 2, 'DisplayName', 'Curva do modelo fracionário Estimado pela Estratégia Computacional 2 \alpha  = 1.10835');
plot(t, obitos_reais, 'ko', 'MarkerSize', 6, 'DisplayName', 'Dados de mortes acumulados de Covid-19 em Bauru'); 
xlabel('Dias (desde 30/03/2020)', 'FontSize', 11);
ylabel('Óbitos acumulados', 'FontSize', 11);
title('Mortes acumuladas de Covid-19 em Bauru', 'FontSize', 14);
ylim([0 30]);  % <-- Aqui você fixa os limites do eixo Y
legend('Location', 'northwest', 'FontSize', 12); grid on;

