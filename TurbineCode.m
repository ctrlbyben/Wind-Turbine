%% Wind Turbine %%
clear; clc; close all;

% Load values in

turbine.profData = readtable("BladeProfile.csv");
turbine.DU91_W2_250 = readtable("DU91-W2-250.csv");
turbine.DU93_W_210 = readtable("DU93-W-210.csv");
turbine.DU96_W_180 = readtable("DU96-W-180.csv");
turbine.DU97_W_300 = readtable("DU97-W-300.csv");


%% Constants 
parameters.rho_air = 1.1;         % Air density [kg/m^3]
parameters.g = 9.81;                % Gravity [m/s^2]
parameters.mu_air = 1.8 * 10^-5;     % [Ns/m^2]

turbine.alpha = 1/3;                %axial induction factor
turbine.nBlades = 3;
turbine.R = max(turbine.profData.DistanceFromCenterOfRotation)/1000; %radius in M
turbine.SweptArea = pi * turbine.R^2;            % Rotor swept area [m^2]
turbine.maxRPM = 15.5; %max turbine speed, RPM
turbine.minRPM = 9.6; %min turbine RPM

%turbine.P = 2.5e6;               % Rated power [W]
%turbine.R = 48;                  % Rotor radius [m]
%turbine.A = pi * (max(R))^2;            % Rotor swept area [m^2]
%turbine.H = 80.4;                % Hub height [m]




%% DELIVERABLE 1

%%Analysis Parameters
nPoints = 101; %number of discrete points for blade analysis
turbine.pitch = 0;   %pitch = 0 deg for first deliverable
turbine.U = 10;      %wind speed = 10m/s for first deliverable
turbine.W = 14*2*pi/60;      %Rad/s, for first deliverable

% Generate Radial Vector (in mm)
r = linspace(min(turbine.profData{:,"DistanceFromCenterOfRotation"}),...
    max(turbine.profData{:,"DistanceFromCenterOfRotation"}), nPoints)/1000;

% CALCULATIONS (Most are just for reference, not necessary for function of
% code

%Calculate AoA, Relative wind Speed, relative wind angle. -----------------ONLY FOR PLOT Remove for Final 
[AoA, Urel, windRelAngle] = calcAoA(r, turbine);

%Calculate Section Twist Angle --------------------------------------------ONLY FOR PLOT Remove for Final 
twist = calcTwist(r, turbine);

%Calculate Cd, Cl ---------------------------------------------------------ONLY FOR PLOT Remove for Final 
[Cd, Cl] = calcCoef(r, turbine, parameters);

%Calculate the Loads on one blade -----------------------------------------ONLY FOR PLOT Remove for Final 
[fTheta, fZ, Lift, Drag] = calcForces(r, turbine, parameters);

%calculate the power [W] and Thrust Load [N] ------------------------------ONLY FOR PLOT Remove for Final 
[Power, Thrust] = calcPower(r, turbine, parameters);

%Calculate Cp and Ct
[Q1.C_p, Q1.C_t] = calcCpCt(r, turbine, parameters);

fprintf("Deliverable 1: \nThe Coefficient of Power is: %.4f \n"  + ...
    "The Coefficient of Thrust is: %.4f \n \n", Q1.C_p, Q1.C_t);

%% DELIVERABLE 2
% Team 24 Specific Parameters:
Q2.U = 8.3; %wind speed, [m/s], for second deliverable
Q2.tipSpdRatio = 7.27; %tip speed ratio for second deliverable
Q2.W = Q2.tipSpdRatio*Q2.U/turbine.R; %new rotational speed for deliverable 2, rad/s

%Update Turbine Struct with this new info
turbine.U = Q2.U;
turbine.W = Q2.W;

%Sweep Through Pitch Angle
Q2.pitchSweep = linspace(-20, 25, 100); %Note: range narrowed down after seeing which pitches cant be evaluated due to lack of AoA-Cd data
Q2.CpSweep = zeros(1, length(Q2.pitchSweep)); %preallocate Cp vector

for ii = 1:length(Q2.pitchSweep)
    turbine.pitch = Q2.pitchSweep(ii);
    [C_p, ~] = calcCpCt(r, turbine, parameters);
    Q2.CpSweep(ii) = C_p;
end

%Plot Cp vs. Pitch Angle
figure(1)
plot(Q2.pitchSweep, Q2.CpSweep)
xlabel("Pitch Angle")
ylabel("C_p")
title("Coefficient of Power vs. Pitch Angle")

Q2.maxCp = max(Q2.CpSweep);
Q2.maxIdx = find(Q2.CpSweep == Q2.maxCp);
Q2.maxCpPitch = Q2.pitchSweep(1,Q2.maxIdx);

fprintf("Deliverable 2: \nThe maximum C_p for a Wind Speed of %.2f [m/s]" + ...
    " and a tip speed ratio of %.2f [-] is %.4f and it occurs at Pitch" + ...
    " = %.2f [deg]. \n \n", Q2.U, Q2.tipSpdRatio, Q2.maxCp, Q2.maxCpPitch);

%% DELIVERABLE 3
%Team 24 Specific Parameters:
Q3.U = 5; %Deliverable 3 wind speed, [m/s]
Q3.pitchSweep = linspace(-25 ,40, 20); %sweep range of pitch
Q3.minTPSPD = turbine.minRPM*2*pi / 60 * turbine.R / Q3.U; %calculate the min tip speed ratio from the min turbine RPM
Q3.maxTPSPD = turbine.maxRPM*2*pi / 60 * turbine.R / Q3.U; %calculate the max tip speed ratio
Q3.tpRatioSweep = linspace(Q3.minTPSPD, Q3.maxTPSPD, 20); %sweep range of tip speed ratio 

%Update Turbine Struct:
turbine.U = Q3.U;

%preallocate array for 2d Cp data
Q3.CpArray = zeros(length(Q3.pitchSweep), length(Q3.tpRatioSweep));

%perform 2d Sweep
for ii = 1:length(Q3.pitchSweep)
    %update pitch angle in turbine struct for this iteration
    turbine.pitch = Q3.pitchSweep(ii);
    for jj = 1:length(Q3.tpRatioSweep)
        %get tip speed ratio for this iteration
        currentTipSpeedRatio = Q3.tpRatioSweep(jj);

        %calculate rotational speed, rad/s
        currentW = currentTipSpeedRatio * Q3.U / turbine.R;

        %set rotational speed in turbine struct for this iteration
        turbine.W = currentW;
        [Cp, ~] = calcCpCt(r, turbine, parameters);

        %store value of Cp in array
        Q3.CpArray(ii,jj) = Cp;
    end
end

figure(2)
surf(Q3.tpRatioSweep, Q3.pitchSweep, Q3.CpArray);
xlabel("Tip Ratio [-]")
ylabel("Pitch Angle [deg]")
zlabel("C_p");

%% Deliverable 4:
%Team 24 Specific Parameters:
Q4.U = 14; %wind speed, m/s for deliverable 4
Q4.W = turbine.maxRPM * 2*pi/60;

%update turbine struct with new data
turbine.U = Q4.U;
turbine.W = Q4.W;

%set up sweep
Q4.pitchsweep = linspace(-20, 25, 100);

Q4.power = zeros([1, length(Q4.pitchsweep)]);%preallocate

for ii = 1:length(Q4.pitchsweep)
    %update pitch value
    turbine.pitch = Q4.pitchsweep(ii);

    %calculate power
    [Q4.power(ii), ~] = calcPower(r, turbine, parameters);
    
end

figure(3)
plot(Q4.pitchsweep, Q4.power);

xlabel("Pitch Angle");
ylabel("Power Output");



%notes: sweep the pitch angle for the given U above, and the max rotational
%speed which is 15.5RPM. Calculate the power for each pitch angle, and find
%the pitch angle that produces rated power. Plot the power vs. wind speed


%% DELIVERABLE 5:
%Tower Analysis: 




%% ----------------Sanity Check Plots-------------------------------------- REMOVE FOR FINAL REPORT
figure(4)
% --- Subplot 1: Aerodynamic Coefficients ---
subplot(3,1,1)
hold on
plot(r, Cl, 'b', 'LineWidth', 1.4)
plot(r, Cd, 'g', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Coefficient Value')
title('Aerodynamic Coefficients Along Blade Span')
legend('C_L','C_D','Location','best')
grid on

% --- Subplot 2: Angles ---
subplot(3,1,2)
hold on
plot(r, windRelAngle, '-b', 'LineWidth', 1.4)
plot(r, twist, '-r', 'LineWidth', 1.4)
plot(r, AoA, '-g', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Angle [deg]')
title('Flow Angles Along Blade Span')
legend('\phi (flow angle)','\theta (twist)','\alpha (AoA)','Location','best')
grid on

% --- Subplot 3: Forces ---
subplot(3,1,3)
hold on
plot(r, Lift, '-b', 'LineWidth', 1.4)
plot(r, Drag, '-r', 'LineWidth', 1.4)
plot(r, fZ, '-g', 'LineWidth', 1.4)
plot(r, fTheta, '-c', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Load [N/m]')
title('Aerodynamic Loads Along Blade Span')
legend('Lift','Drag','F_z','F_\Theta')
grid on



%% -------------FUNCTIONS:------------------ %%

%% Calculate Section Twist for every r
function [sectionTwist] = calcTwist(r, turbine)
    r = r*1000; %convert r to mm
    BladeData = turbine.profData;
    sectionTwist = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.BladeTwist, r);
end

%% Calculate Chord Length for every r
function [chordLength] = calcChord(r, turbine)
    r = r*1000; %convert r to mm
    BladeData = turbine.profData;
    chordLength = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.ChordLength, r);
    chordLength = chordLength / 1000; %return chord in meters
end

%% Calculate Angular Induction Factor for every r
function alphaPrime = calcAPrime(r, turbine)

    W = turbine.W; %turbine rotational speed, rad/s
    lambdaR = (W .*r)/turbine.U; %get local speed ratio
    alpha = turbine.alpha; %get axial induction factor

    alphaPrime = -1/2 + 1/2.*sqrt(1+4./((lambdaR).^2).*alpha.*(1-alpha)); %calculate angular induction factor
  
end

%% Calculate Angle of Attack for every r
function [AoA,Urel,windRelAngle]  = calcAoA(r,turbine)
%CALCAoA calculates the angle of attack of an airfoil

    W = turbine.W; %turbine rotational speed in rad/s
    U = turbine.U; %wind speed in m/s

    a = turbine.alpha; %get a value
    a_prime = calcAPrime(r,turbine); %get a prime value vector

    sectionPitch = calcTwist(r, turbine) + turbine.pitch; %calculate section pitch

    %calculate the angle of relative wind and the vertical
    tangentComponent = r.*W.*(1+a_prime);
    axialComponent = (U.*(1-a));

    windRelAngle = atand(axialComponent ./ tangentComponent);

    %calculate the angle between the 
    AoA = windRelAngle-sectionPitch;
    
    %calculate magnitude of Urel
    Urel = sqrt(axialComponent.^2 + tangentComponent.^2);
end

%% Calculate airfoil type for every r
function [airfoil] = getAirfoil(r)
    r = r*1000; %convert to mm
    airfoil = strings(1,length(r)); %preallocate array

    for ii = 1:length(r) %filter r value to get airfoil type
        if r(ii) < 11625
            airfoil(ii) = "circle";
        elseif r(ii)<23250
            airfoil(ii) = "DU 97-W-300";
        elseif r(ii)<34875
            airfoil(ii) = "DU 91-W2-250";
        elseif r(ii)<44175
            airfoil(ii) = "DU 93-W-210";
        elseif r(ii)<=48000
            airfoil(ii) = "DU 96-W-180";
        end
    end
end

%% Calculate Cd and Cl for every r
function [Cd, Cl] = calcCoef(r, turbine, parameters)
    
    %first get the airfoil type for each r (use if-then statements on r and blade profile data)
    airfoils = getAirfoil(r);

    %then get the AoA of each R using function
    [AoA,Urel,~] = calcAoA(r, turbine);

    %get chord lengths (used for Cd of circle)
    chords = calcChord(r, turbine);

    %using AoA and airfoil type for each r, get coefficients for each R
    Cl = zeros(1,length(r)); %preallocate
    Cd = zeros(1,length(r)); %preallocate

    %Get Cl and Cd value depending on airfoil type and interpolating for AoA
    for ii = 1:length(r)

        localAirfoil = airfoils(ii);
        localAoA = AoA(ii);

        if localAirfoil == "circle"
            Cl(ii) = 0;
            %calculate reynolds number of section
            Re = parameters.rho_air .* Urel(ii) .* chords(ii) / parameters.mu_air;
            
            %calculate CD using external function
            Cd(ii) = cylinderCD(Re); 

        elseif localAirfoil == "DU 97-W-300"
            coefTable = turbine.DU97_W_300;
            Cl(ii) = interp1(coefTable.AoA, coefTable.CL, localAoA);
            Cd(ii) = interp1(coefTable.AoA, coefTable.CD, localAoA);
        elseif localAirfoil == "DU 91-W2-250"
            coefTable = turbine.DU91_W2_250;
            Cl(ii) = interp1(coefTable.AoA, coefTable.CL, localAoA);
            Cd(ii) = interp1(coefTable.AoA, coefTable.CD, localAoA);
        elseif localAirfoil == "DU 93-W-210"
            coefTable = turbine.DU93_W_210;
            Cl(ii) = interp1(coefTable.AoA, coefTable.CL, localAoA);
            Cd(ii) = interp1(coefTable.AoA, coefTable.CD, localAoA);
        elseif localAirfoil == "DU 96-W-180"
            coefTable = turbine.DU96_W_180;
            Cl(ii) = interp1(coefTable.AoA, coefTable.CL, localAoA);
            Cd(ii) = interp1(coefTable.AoA, coefTable.CD, localAoA);
        end
    end

end

%% Calculate Drag and Lift Forces
function [fTheta, fZ, Lift, Drag] = calcForces(r, turbine, parameters)
    [~,Urel,windRelAngle]  = calcAoA(r,turbine); % get relative wind speed and angle (deg)
    chordLength = calcChord(r,turbine); %get chord length for every r
    [Cd, Cl] = calcCoef(r, turbine, parameters); %get coefficients
    
    
    Lift = 0.5 .* parameters.rho_air .* Urel.^2 .* Cl .* chordLength;
    Drag = 0.5 * parameters.rho_air .* Urel.^2 .* Cd .* chordLength;

    fTheta = Lift .* sind(windRelAngle) - Drag .* cosd(windRelAngle);
    fZ = Lift .*cosd(windRelAngle) + Drag.* sind(windRelAngle);

end

%% Calculate Power the turbine

function [Power, Thrust] = calcPower(r, turbine, parameters)
    [fTheta, fZ, ~, ~] = calcForces(r, turbine, parameters);
    
    d_Torque = fTheta .* r;
    Torque = trapz(r, d_Torque);
    Power = turbine.nBlades* Torque * turbine.W; 

    Thrust = turbine.nBlades* trapz(r, fZ);

    
end

%% Calculate Cp and Ct
function [C_p, C_t] = calcCpCt(r, turbine, parameters)
    [Power, Thrust] = calcPower(r, turbine, parameters);

    %calculate coefficients
    C_p = Power / (1/2 * parameters.rho_air * turbine.U^3 ...
        * turbine.SweptArea);
    C_t = Thrust / (1/2 * parameters.rho_air * turbine.U^2 ...
        * turbine.SweptArea);
end
