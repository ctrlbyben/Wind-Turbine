%% Wind Turbine %%
clear; clc; close all;

% Load values in

BladePitch = 0; %(0 in Deliverable 1)
Blade = readtable('Turbine Data/BladeProfile.csv');

% Remove hub/circle sections
Blade = Blade(~strcmpi(Blade.Airfoil, 'circle'), :);

% Radial Vector
r = Blade.DistanceFromCenterOfRotation;   % vector of all radii [mm]
r = r / 1000; % Converting from [mm] to [m]

theta = Blade.BladeTwist; % Blade Twist from Blade Profile 
beta = BladePitch; % Blade Pitch (0 in Deliverable 1)
nSections = height(Blade);


%% Constants 
V_infin = 10; %[m/s]
W = 14 * (2*pi / 60); % [Rad/s]
rho_air = 1.1;         % Air density [kg/m^3]
g = 9.81;                % Gravity [m/s^2]
mu_air = 1.8 * 10^-5;     % [Ns/m^2]
P = 2.5e6;               % Rated power [W]
R = 48;                  % Rotor radius [m]
A = pi * (max(r))^2;            % Rotor swept area [m^2]
H = 80.4;                % Hub height [m]
NumberOfBlades = 3;

% Preallocating angles and aerodynamic coefficient vectors
phi = zeros(nSections,1);
alpha = zeros(nSections,1);
V_rel = zeros(nSections,1);
CL = zeros(nSections,1);
CD = zeros(nSections,1);
CM = zeros(nSections,1);

for i = 1:nSections
    % Local Relative Velocity
    V_rel(i) = sqrt(V_infin^2 + (W * r(i))^2);
    
    % Angle of Relative Wind [rad]
    phi(i) = atan(V_infin / (W * r(i)));

    % Angle of Attack
    alpha(i) = rad2deg(phi(i)) - (theta(i)+ beta);
end

for i = 1:nSections
    [CL(i), CD(i), CM(i), sectionInfo] = getBladeSectionCharacteristics(r(i), alpha(i));
end

%Plot CL and CD vs Radius for quality check
%% --- Subplot 1: Aerodynamic Coefficients ---
subplot(2,1,1)
hold on
plot(r, CL, 'b', 'LineWidth', 1.4)
plot(r, CD, 'g', 'LineWidth', 1.4)
plot(r, CM, 'r', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Coefficient Value')
title('Aerodynamic Coefficients Along Blade Span')
legend('C_L','C_D','C_M','Location','best')
grid on

%% --- Subplot 2: Angles ---
subplot(2,1,2)
hold on
plot(r, rad2deg(phi), '-b', 'LineWidth', 1.4)
plot(r, theta, '-r', 'LineWidth', 1.4)
plot(r, alpha, '-g', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Angle [deg]')
title('Flow Angles Along Blade Span')
legend('\phi (flow angle)','\theta (twist)','\alpha (AoA)','Location','best')
grid on
%% Compute Lift and Drag forces on the blade
% Tangential and Axial components

% Preallocate D and L Vectors

D = zeros(nSections,1);
L = zeros(nSections,1);

ChordLength = Blade.ChordLength;

% Compute Lift and Drag forces on the blade
for i = 1:nSections
    L(i) = 0.5 * rho_air * V_rel(i)^2 * CL(i) * ChordLength(i);
    D(i) = 0.5 * rho_air * V_rel(i)^2 * CD(i) * ChordLength(i);
end

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

figure;
subplot(2,1,1);
plot(r .*10^-3, Blade.ChordLength,'LineWidth',1.5);
xlabel('Distance from Center [m]');
ylabel('Chord Length [m]');
title('Chord Length vs Center of Rotation');

subplot(2,1,2);
plot(Blade.DistanceFromCenterOfRotation, Blade.BladeTwist,'LineWidth',1.5);
xlabel('Distance from Center [m]');
ylabel('Blade Twist [degrees]');
title('Blade Twist vs Center of Rotation');


%% TO BE USED LATER %% 
% V_cut_in  = 4;           % Cut-in wind speed [m/s]
% V_rated   = 11;          % Rated wind speed [m/s]
% V_cut_out = 25;          % Cut-out wind speed [m/s]
% omega_min = 9.6 * 2*pi/60;  % Minimum rotor speed [rad/s]
% omega_max = 15.5 * 2*pi/60; % Maximum rotor speed [rad/s]

% % Tip-speed ratio [ ]
% lambda = (omega * R) / V;

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
