function [Pxx_slope_LS,Pxx_line_LS] = spectra_linefit2(W_LS,Pxx_LS_mean,mindX,maxdX)

ft = fittype( 'poly1' );
% mindX = 2e-5;maxdX = 3e-4;

Jl = find(mindX < W_LS & maxdX > W_LS);
[model, ~] = fit(log(W_LS(Jl)),log(Pxx_LS_mean(Jl)), ft );
coeff_vals = coeffvalues(model);
Pxx_slope_LS = coeff_vals(1); 
Pxx_line_LS = model(log(W_LS(Jl)));



