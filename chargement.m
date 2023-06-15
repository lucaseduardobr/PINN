function load_f = chargement(X,x_F_min,x_F_max,F,w_max,numInternalCollocationPoints)

porte=@(x) x<x_F_max & x>x_F_min ;

load_f = porte(X); % on divise par w_max à cause de la normalisation de W

l = x_F_max-x_F_min;

%     load_f = F/l*load_f/w_max; %

load_f = F/l*load_f;

end