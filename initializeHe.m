function parameter = initializeHe(sz,numIn,className)
%rng(1234);
arguments
    sz
    numIn
    className = 'single'
end

parameter = sqrt(2/numIn) * randn(sz,className);
parameter = dlarray(parameter);

end