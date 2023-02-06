function time = benchmark_efms(S)
    rev = zeros(1,size(S,2));
    tic
    calculate_flux_modes(S,rev);
    time = toc;
end