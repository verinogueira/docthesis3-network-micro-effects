%% FUNCTION PERFORM APPROXIMATION EXERCISE FOR LEONTIEF INVERSE
% n - limit of the order of approximation
% N - dimension of Gamma
% Gamma - the Gamma matrix
% shock - the industry-specific shock vector


function ApEf = AppEffect(n, N, Gamma, shock)
    if n==0;
    ApEf = eye(N)*shock;
    end
     
        part=0;
    for i = 1:n
        part = part + Gamma'^i;
    end
     ApEf = (eye(N) + part) * shock;
end