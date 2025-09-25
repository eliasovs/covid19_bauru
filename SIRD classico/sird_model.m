function dYdt = sird_model(t, Y, beta, lambda, gamma)
    S = Y(1);
    I = Y(2);
    R = Y(3);
    D = Y(4);

    dSdt = -beta * S * I;
    dIdt = beta * S * I - lambda * I - gamma * I;
    dRdt = lambda * I;
    dDdt = gamma * I;

    dYdt = [dSdt; dIdt; dRdt; dDdt];
end
