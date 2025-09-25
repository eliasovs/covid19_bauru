
function res = func_residuo_frac_multi_alpha(param, t0, tf, h, y0, casos_reais, obitos_reais)
    beta = param(1); lambda = param(2); gamma = param(3);
    alpha1 = param(4); alpha2 = param(5); alpha3 = param(6); alpha4 = param(7);

    f_fun = @(t, y, p) [
        -(p(1)^alpha1) * y(1) * y(2);
         (p(1)^alpha2) * y(1) * y(2) - (p(2)^alpha2) * y(2) - (p(3)^alpha2) * y(2);
         (p(2)^alpha4) * y(2);
         (p(3)^alpha3) * y(2)
    ];

    S0 = y0(1); I0 = y0(2);
    dy0 = [
        -beta * S0 * I0;
         beta * S0 * I0 - lambda * I0 - gamma * I0;
         lambda * I0;
         gamma * I0
    ];
    y_inicial = [y0, dy0];

    [~, y] = fde_pi12_pc([alpha1, alpha2, alpha3, alpha4], f_fun, t0, tf, y_inicial, h, [beta, lambda, gamma]);

    I = y(2,:); R = y(3,:); D = y(4,:);
    casos_modelo = I + R + D;
    obitos_modelo = D;

    res = [casos_modelo(:) - casos_reais(:); obitos_modelo(:) - obitos_reais(:)];
end
