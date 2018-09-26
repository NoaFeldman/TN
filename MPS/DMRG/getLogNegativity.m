function epsilon = getLogNegativity(rhoT2)
    epsilon = log(sum(abs(getSpec(rhoT2))));
end
