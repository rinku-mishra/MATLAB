clc; clearvars;
[ks1,X11]=surface_plot_wpiA_ratio;
[ks2,X12]=surface_plot_wpi_ratio;
    figure(1)
    subplot(211)
    plot(ks1,X11,'linewidth',2);
    ylabel('\omega/\omega_{piA}');
    xlabel('\kappa\lambda_{eA}');
    grid on
   
    subplot(212)
    plot(ks2,X12,'linewidth',2);
    ylabel('\omega/\omega_{piB}');
    xlabel('\kappa\lambda_{eA}')
    grid on
   
   
%     legend
   
%     ylabel('\omega/\omega_{peB}');
%     xlabel('\kappa\lambda_{eB}');
    grid on
