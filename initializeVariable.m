function parameter = initializeVariable(sz,valeur_variable, className)

arguments
    sz
    valeur_variable
    className = 'single'
end

parameter = valeur_variable*ones(sz,className);
parameter = dlarray(parameter);

end