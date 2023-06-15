function weights = initializeGlorot(sz,numOut,numIn)
rng(1234);
Z = 2*rand(sz,'single') - 1;
bound = sqrt(6 / (numIn + numOut));

weights = bound * Z;
weights = dlarray(weights);

end

