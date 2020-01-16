%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 
%%%%% Reproduction Functions:
%%%%% 
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Beverton-Holt:
% Growth parameters:
  R_0   = 3;             % intrinsic growth rate
  K = 1.0;                 % carrying capacity
% Plot variables:
  N_plot  = [0.0 : 0.05 : K];
  FN_plot = (R_0 * N_plot) ./ (1 + ((R_0 - 1)/K)*N_plot);
%  cstar = sqrt(2*sigma2*log(R_0));

%%%%% Ricker:
% Growth parameters:
 R1 = 1.5;
 K1 = 1.0;
%%% Plot variables:
 N_plot  = [0.0 : 0.05 : K1];
 FN1_plot = N_plot .* exp(R1 * (1 - (N_plot/K1)));
 
 %%%%% Sigmoid Beverton-Holt;
% Growth parameters:
  R2     = 5;
  K2     = 1.0;
  delta = 2;
% % %%% Plot variables:
  N_plot  = [0.0 : 0.05 : K2];
  FN2_plot = (R2 * N_plot.^(delta)) ./ (1 + (K2^(delta-1)*R2-1)/(K2^(delta))*N_plot.^(delta));

%%%%% Plot the growth function:
% subplot(2,1,1);
plot( N_plot, FN_plot, 'k--');
hold on
plot( N_plot, FN1_plot, 'k:');
hold on
plot( N_plot, FN2_plot, 'k-.');
hold on
plot( N_plot, N_plot, 'k-');%,'HandleVisibility','off');
hold off
xlabel( '$y$');
ylabel( '$f(y)$');
title(  'Growth functions');
legend('Beverton-Holt', 'Ricker', 'Sigmoid Beverton-Holt', 'Reference Line')
