
% %% Wind Turbine %%
% 
% clear; 
% clc; 
% close all;
% 
% % Constants
% rho_air = 1.225;         % Air density [kg/m^3]
% g = 9.81;                % Gravity [m/s^2]
% 
% P = 2.5e6;               % Rated power [W]
R = 48;                  % Rotor radius [m]
% A = pi * R^2;            % Rotor swept area [m^2]
% H = 80.4;                % Hub height [m]
% 
% V_cut_in  = 4;           % Cut-in wind speed [m/s]
% V_rated   = 11;          % Rated wind speed [m/s]
% V_cut_out = 25;          % Cut-out wind speed [m/s]
% omega_min = 9.6 * 2*pi/60;  % Minimum rotor speed [rad/s]
% omega_max = 15.5 * 2*pi/60; % Maximum rotor speed [rad/s]
% V_air = ;                   % Wind Velocity
% C_p = ; % Coefficient of Power [ ] Must be <.5926 - Typically .3-.5
% a = ; % Axial Induction Factor [ ] - Typically .2-.4
% 
% % Tip-speed ratio [ ]
% lambda = (omega * R) / V;
% 
% % Power extracted from the wind [W]
% W_dot = 0.5 * rho * A * V_;
% 
% % Useful mechanical power from turbine [W]
% P_turbine = P_wind * C_p;
% 
% % Thrust force on rotor [N]
% Thrust = 0.5 * rho * A * V^2 * CT;
% 
% % Torque on rotor [Nm]
% Torque = P_turbine / omega;
% 
% % % Deliverable 1 % %
V = 10; % [m/s]
V_rotation = 14; % [RPM]
beta = 0; % [Degrees]
omega = ((V_rotation*2*pi)/60); % rad/s
lambda_1 = omega*R/V;

% % Sweep Through a Values in C_p and C_t
a = linspace(0, 0.5, 100);
C_p = 4.*a.*(1 - a).^2;
C_t = 4.*a.*(1 - a);

figure;
subplot(2,1,1)
plot(a, C_p, 'LineWidth',1.5)
xlabel('Axial Induction Factor a')
ylabel('C_p')
title('Power Coefficient vs a')
grid on

subplot(2,1,2)
plot(a, C_t, 'LineWidth',1.5)
xlabel('Axial Induction Factor a')
ylabel('C_t')
title('Thrust Coefficient vs a')
grid on
