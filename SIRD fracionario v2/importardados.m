% === IMPORTAÇÃO DA BASE DE DADOS COVID-19 DE BAURU ===

% Nome do arquivo
arquivo = 'dadoscovidBauru.xlsx';

% Leitura completa da tabela
dados = readtable(arquivo);

% Seleção das colunas relevantes
dias              = dados.Dia;
casos_acumulados  = dados.Casos_Acumulados;
obitos_acumulados = dados.Obitos_Acumulados;
recuperados       = dados.Recuperados;

% Visualização de verificação
disp('Primeiras linhas dos dados importados:');
T = table(dias, casos_acumulados, obitos_acumulados, recuperados);
disp(head(T, 10));  % Exibe as 10 primeiras linhas

% === GRÁFICO: CASOS ACUMULADOS AO LONGO DO TEMPO ===

figure;
plot(dias, casos_acumulados, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Dia');
ylabel('Casos Acumulados');
title('Evolução dos Casos Acumulados de COVID-19 em Bauru');
grid on;
xlim([min(dias) max(dias)]);
ylim([0 max(casos_acumulados)*1.1]);

% === GRÁFICO: CASOS ACUMULADOS - APENAS PONTOS ===

figure;
plot(dias, casos_acumulados, 'o', 'MarkerSize', 6, 'LineWidth', 1.5);
xlabel('Dia');
ylabel('Casos Acumulados');
title('Evolução dos Casos Acumulados de COVID-19 em Bauru');
grid on;
xlim([min(dias) max(dias)]);
ylim([0 max(casos_acumulados)*1.1]);
