%% Wind Turbine %%
clear; clc; close all;

% Load values in

BladePitch = 0; %(0 in Deliverable 1)
Blade = readtable('Turbine Data/BladeProfile.csv');

% Remove hub/circle sections - there not contributing to the power
Blade = Blade(~strcmpi(Blade.Airfoil, 'circle'), :);

% Radial Vector
r = Blade.DistanceFromCenterOfRotation;   % vector of all radii [mm]
r = r / 1000; % Converting from [mm] to [m]

theta = Blade.BladeTwist; % Blade Twist from Blade Profile 
beta = BladePitch; % Blade Pitch (0 in Deliverable 1)
nSections = height(Blade);


%% Constants 
V_infin = 10;            % [m/s]
W = 14 * (2*pi / 60);    % [Rad/s]
rho_air = 1.1;           % Air density [kg/m^3]
g = 9.81;                % Gravity [m/s^2]
mu_air = 1.8 * 10^-5;    % [Ns/m^2]
P = 2.5e6;               % Rated power [W]
R = 48;                  % Rotor radius [m]
A = pi * (max(r))^2;     % Rotor swept area [m^2]
H = 80.4;                % Hub height [m]
NumberOfBlades = 3;

% Preallocating angles and aerodynamic coefficient vectors
phi = zeros(nSections,1); % Angle of Relative Wind
alpha = zeros(nSections,1); % Angle of Attack
V_rel = zeros(nSections,1); % Relative Velocity
CL = zeros(nSections,1); % Coefficient of Lift
CD = zeros(nSections,1); % Coefficient of Drag
CM = zeros(nSections,1); % Coefficient of Moment

for i = 1:nSections
    % Local Relative Velocity
    V_rel(i) = sqrt(V_infin^2 + (W * r(i))^2);
    
    % Local Angle of Relative Wind [rad]
    phi(i) = atan(V_infin / (W * r(i)));

    % Local Angle of Attack
    alpha(i) = rad2deg(phi(i)) - (theta(i)+ beta);
end

for i = 1:nSections % Finds each coefficient for each section of blade and AoA
    [CL(i), CD(i), CM(i), sectionInfo] = getBladeSectionCharacteristics(r(i), alpha(i));
end

% Plot CL and CD vs Radius for quality check
% Subplot 1: Aerodynamic Coefficients
subplot(2,1,1)
hold on
plot(r, CL, 'b', 'LineWidth', 1.5)
plot(r, CD, 'g', 'LineWidth', 1.5)
% plot(r, CM, 'r', 'LineWidth', 1.5)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Coefficient Value')
title('Aerodynamic Coefficients Along Blade Span')
legend('C_L','C_D','Location','best')
grid on

% Subplot 2: Angles
subplot(2,1,2)
hold on
plot(r, rad2deg(phi), '-b', 'LineWidth', 1.5)
plot(r, theta, '-r', 'LineWidth', 1.5)
plot(r, alpha, '-g', 'LineWidth', 1.5)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Angle [deg]')
title('Flow Angles Along Blade Span')
legend('\phi (Relative Wind Angle)','\theta (Blade Twist)','\alpha (AoA)','Location','best')
grid on

%% Compute dift and drag forces on the blade %%

% Preallocate dD and dL Vectors - Chord Area Equation requires chord length
% and dr

dr = diff(r); % Takes the difference between each element of r
dr = [dr; dr(end)]; % Makes diff(r) the same length as the original r

% Preallocating dD and dL
dD = zeros(nSections,1);
dL = zeros(nSections,1);

ChordLength = Blade.ChordLength / 1000; % [mm] -> [m]

% Compute Lift and Drag forces on each blade section (dr)
for i = 1:nSections
    dL(i) = 0.5 * rho_air * V_rel(i)^2 * CL(i) * ChordLength(i) * dr(i); % Lift on each blade section
    dD(i) = 0.5 * rho_air * V_rel(i)^2 * CD(i) * ChordLength(i) * dr(i); % Drag on each blade section 
end
% Tangential and Axial components (torque and thrust)
dF_Thrust = (dL .* cos(phi) + dD .* sin(phi));  % Thrust force = number of blades * horizonal comp of lift force and veritical of drag 
dF_Torque = (dL .* sin(phi) - dD .* cos(phi));
 
%Plot thrust and torque 

TotalThrust = sum(dF_Thrust);
TotalTorque = sum(r .* dF_Torque); % Trapz Instead of sum for both

C_P_ = ((TotalTorque * W) / (.5*(rho_air*A*V_infin^3)));
C_T_ = (TotalThrust / (.5*rho_air*A*V_infin^2));

% Axial induction factor from thrust coefficient
a = (1 - sqrt(1 - C_T_)) / 2;

fprintf('Axial induction factor a = %.4f\n', a);
CP_check = 4 * a * (1 - a)^2;
fprintf('Predicted CP (theoretical) = %.4f\n', CP_check);


% % Sweep Through a Values in C_p and C_t
a = linspace(0, 0.5, 100);
C_p = 4.*a.*(1 - a).^2;
C_t = 4.*a.*(1 - a);

% W_given = ; % Will be group specific
% V_given = ; % Will be group specific

% % Tip Speed Ratio
% lambda_given = ;% Will be group specific
% W_given = ;
% (W_given * R) / V_given; 
% 
V_infin = 8; % Examples until we get the real one
lambda = 6; % Example until we get the real one 
W = lambda * V_infin / R;   % [rad/s]

% Sweeping through beta values
beta_array = -5:.5:15;        % pitch angles to test [deg]
CP_array_beta = zeros(size(beta_array));
CT_array_beta = zeros(size(beta_array));

for b = 1:length(beta_array)
    beta = beta_array(b); % Currgent pitch angle
    % Recompute angles & section coefficients for this beta
    for i = 1:nSections
        phi(i) = atan(V_infin / (W * r(i)));            % [rad]
        alpha(i) = rad2deg(phi(i)) - (theta(i) + beta); % [deg]
        [CL(i), CD(i), ~] = getBladeSectionCharacteristics(r(i), alpha(i));
    end
     dL = 0.5 * rho_air .* V_rel.^2 .* CL .* ChordLength .* dr;
     dD = 0.5 * rho_air .* V_rel.^2 .* CD .* ChordLength .* dr;

    % Tangential and Axial Components
    dF_Thrust = (dL .* cos(phi) + dD .* sin(phi));  
    dF_Torque = (dL .* sin(phi) - dD .* cos(phi));

    % Integrating to find total thrust and torque
    TotalTorque = sum(r .* dF_Torque);
    TotalThrust = sum(dF_Thrust);


    P_total = W * TotalTorque;    % Total Power = Omega x TotalTorque
    A = pi * (R)^2; % Swept Area

    CP_array_beta(b) = P_total / (0.5 * rho_air * A * V_infin^3);
    CT_array_beta(b) = TotalThrust / (0.5 * rho_air * A * V_infin^2);

end

% Plotting the sweep of beta values vs C_P
figure;
plot(beta_array,CP_array_beta,'LineWidth',1.5);
hold on
xlabel('Pitch Angle [degree]');
ylabel('Coefficient of Power');
title('Pitch Angle vs Coefficient of Power');
[maxY, maxIndex] = max(CP_array_beta);

% Mark the maximum value with a red dot
plot(beta_array(maxIndex), maxY, 'r', 'MarkerSize', 10, 'DisplayName', 'Max Cp');
% Label the maximum point with x and y values
text(beta_array(maxIndex), maxY, sprintf(' (x: %d, y: %.4f)', beta_array(maxIndex), maxY), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red');


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
plot(r, Blade.ChordLength,'LineWidth',1.5);
xlabel('Distance from Center [m]');
ylabel('Chord Length [m]');
title('Chord Length vs Center of Rotation');

subplot(2,1,2);
plot(r, Blade.BladeTwist,'LineWidth',1.5);
xlabel('Distance from Center [m]');
ylabel('Blade Twist [degrees]');
title('Blade Twist vs Center of Rotation'); % Chord length gets smaller radially - Because the fin moves faster radially


%% TO BE USED LATER %% 
% V_cut_in  = 4;           % Cut-in wind speed [m/s]
% V_rated   = 11;          % Rated wind speed [m/s]
% V_cut_out = 25;          % Cut-out wind speed [m/s]
% omega_min = 9.6 * 2*pi/60;  % Minimum rotor speed [rad/s]
% omega_max = 15.5 * 2*pi/60; % Maximum rotor speed [rad/s]



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
