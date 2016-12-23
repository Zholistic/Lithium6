close all;
dipole_firstmodes = [22.24, 22.4, 22.18, 22.98, 22.9, 24.55, 24.82];
dipole_firstmodes_error = [0.045, 0.07, 0.03, 0.04, 0.16, 0.09, 0.04];

mode_fields = [680, 690, 700, 735, 750, 832, 860];

breathing_modes_pca = [45.2, 49.91, NaN, 46.78, NaN, 49.14, 50.03];
breathing_modes_pca_error = [0.27, 1, NaN, 0.19, NaN, 0.21, 0.63];

breathing_modes_radialfit = [NaN, NaN, 46.36, NaN, 47.32, NaN, NaN];
breathing_modes_radialfit_error = [NaN, NaN, 0.17, NaN, 0.32, NaN, NaN];

%Error calculation:
deltaOmega_b_on_omega_o = (breathing_modes_pca - 2*dipole_firstmodes)./dipole_firstmodes;
deltaOmega_b_on_omega_o_radial = (breathing_modes_radialfit - 2*dipole_firstmodes)./dipole_firstmodes;

numerator_error = breathing_modes_pca_error + dipole_firstmodes_error;
numerator_error_radial = breathing_modes_radialfit_error + dipole_firstmodes_error;
deltaOmega_b_on_omega_o_error = deltaOmega_b_on_omega_o.*((numerator_error./(breathing_modes_pca - 2*dipole_firstmodes)) + dipole_firstmodes_error./dipole_firstmodes);
deltaOmega_b_on_omega_o_error_radial = deltaOmega_b_on_omega_o_radial.*((numerator_error_radial./(breathing_modes_radialfit - 2*dipole_firstmodes)) + dipole_firstmodes_error./dipole_firstmodes);
%---------------------%

figure(1);
errorbar(mode_fields,dipole_firstmodes,dipole_firstmodes_error,'.');
grid on;
title(['Dipole frequency (first PCA mode) vs Field (G)']); 

figure(2);
errorbar(mode_fields,breathing_modes_pca,breathing_modes_pca_error,'.');
hold on;
errorbar(mode_fields,breathing_modes_radialfit,breathing_modes_radialfit_error,'.');
grid on;
title(['Breathing frequency (PCA mode where possible else radial width) vs Field (G)']);

figure(3);
errorbar(mode_fields,(breathing_modes_pca - 2*dipole_firstmodes)./dipole_firstmodes,deltaOmega_b_on_omega_o_error,'.');
hold on;
errorbar(mode_fields,(breathing_modes_radialfit - 2*dipole_firstmodes)./dipole_firstmodes,deltaOmega_b_on_omega_o_error_radial,'.');
grid on;
axis([650 900 -0.05 0.3]);
title(['\delta \omega_B on \omega_o vs Field (G)']);
