%% === MÉTRICAS DE AVALIAÇÃO PARA O MODELO SIRD FRACIONÁRIO ===

% Vetores reais e simulados (62 dias usados)
casos_reais_frac   = casos_acumulados(1:62);
casos_estimados_frac = casos_modelo(1:62);

obitos_reais_frac  = obitos_acumulados(1:62);
obitos_estimados_frac = obitos_modelo(1:62);

% Garante formato coluna
casos_estimados_frac  = casos_estimados_frac(:);
casos_reais_frac      = casos_reais_frac(:);
obitos_estimados_frac = obitos_estimados_frac(:);
obitos_reais_frac     = obitos_reais_frac(:);

%% === EQM e REQM ===
eqm_casos_frac  = mean((casos_reais_frac - casos_estimados_frac).^2);
eqm_obitos_frac = mean((obitos_reais_frac - obitos_estimados_frac).^2);

reqm_casos_frac  = sqrt(eqm_casos_frac);
reqm_obitos_frac = sqrt(eqm_obitos_frac);

fprintf('Modelo fracionário com alpha - 1.10835, 62 dias');
fprintf('\n=== EQM e REQM (Modelo Fracionário) ===\n');
fprintf('EQM - Casos acumulados  = %.4f\n', eqm_casos_frac);
fprintf('REQM - Casos acumulados = %.4f\n', reqm_casos_frac);
fprintf('EQM - Óbitos acumulados  = %.4f\n', eqm_obitos_frac);
fprintf('REQM - Óbitos acumulados = %.4f\n', reqm_obitos_frac);

%% === MAPE ===

idx_valid_casos_frac  = casos_reais_frac > 0;
idx_valid_obitos_frac = obitos_reais_frac > 0;

mape_casos_frac = mean( ...
    abs(casos_reais_frac(idx_valid_casos_frac) - casos_estimados_frac(idx_valid_casos_frac)) ...
    ./ casos_reais_frac(idx_valid_casos_frac) );

mape_obitos_frac = mean( ...
    abs(obitos_reais_frac(idx_valid_obitos_frac) - obitos_estimados_frac(idx_valid_obitos_frac)) ...
    ./ obitos_reais_frac(idx_valid_obitos_frac) );

fprintf('\n=== MAPE (Modelo Fracionário) ===\n');
fprintf('MAPE - Casos acumulados  = %.4f\n', mape_casos_frac);
fprintf('MAPE - Óbitos acumulados = %.4f\n', mape_obitos_frac);

%% === ICC ===

M_frac_casos = [casos_reais_frac(:), casos_estimados_frac(:)];
[r_frac, LB_frac, UB_frac, F_frac, df1_frac, df2_frac, p_frac] = ICC(M_frac_casos, 'A-1', 0.05, 0);

fprintf('\n=== ICC para Casos Acumulados (Fracionário) ===\n');
fprintf('ICC       = %.4f\n', r_frac);
fprintf('IC 95%%    = [%.4f, %.4f]\n', LB_frac, UB_frac);
fprintf('F         = %.4f\n', F_frac);
fprintf('GL df1    = %d\n', df1_frac);
fprintf('GL df2    = %.2f\n', df2_frac);
fprintf('p-valor   = %.4f\n', p_frac);

M_frac_obitos = [obitos_reais_frac(:), obitos_estimados_frac(:)];
[r_frac_ob, LB_frac_ob, UB_frac_ob, F_frac_ob, df1_frac_ob, df2_frac_ob, p_frac_ob] = ICC(M_frac_obitos, 'A-1', 0.05, 0);

fprintf('\n=== ICC para Óbitos Acumulados (Fracionário) ===\n');
fprintf('ICC       = %.4f\n', r_frac_ob);
fprintf('IC 95%%    = [%.4f, %.4f]\n', LB_frac_ob, UB_frac_ob);
fprintf('F         = %.4f\n', F_frac_ob);
fprintf('GL df1    = %d\n', df1_frac_ob);
fprintf('GL df2    = %.2f\n', df2_frac_ob);
fprintf('p-valor   = %.4f\n', p_frac_ob);

%%
%% === EQM e REQM TOTAIS (Casos + Óbitos) ===
residuos_totais = [casos_reais_frac - casos_estimados_frac; obitos_reais_frac - obitos_estimados_frac];
eqm_total_frac = mean(residuos_totais.^2);
reqm_total_frac = sqrt(eqm_total_frac);

fprintf('\n=== EQM e REQM TOTAIS (Casos + Óbitos) ===\n');
fprintf('EQM Total  = %.4f\n', eqm_total_frac);
fprintf('REQM Total = %.4f\n', reqm_total_frac);
