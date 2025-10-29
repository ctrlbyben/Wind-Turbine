%% ME5053 PROJECT 2
%Wind Turbine Analysis
%Ben Lahyani, Simon Granberg, Austin Gage

clear; clc; close all;

%% LOAD DATA

turbine.profileData = readtable("BladeProfile.csv");
turbine.DU91_W2_250 = readtable("DU91-W2-250.csv");
turbine.DU93_W_210 = readtable("DU93-W-210.csv");
turbine.DU96_W_180 = readtable("DU96-W-180.csv");
turbine.DU97_W_300 = readtable("DU97-W-300.csv");
turbine.towerSpecs = readtable("towerSpecs.csv", ...
    VariableNamingRule="preserve"); 

%% PARAMETERS
param.rho_air = 1.1;                 % Air density [kg/m^3]
param.g = 9.81;                      % Gravity [m/s^2]
param.mu_air = 1.8 * 10^-5;          % Air Dynamic Viscosity [Ns/m^2]
param.rho_steel = 7850;              % steel density [kg/m^3]
param.TS_steel = 450;                % tensile strength of steel [MPa]
param.sigmaY_steel = 345;            % yield strength of steel [MPa]
param.E_steel = 200;                 % youngs modulus of steel  [GPa]
param.Cl = 1;                        % fatigue loading factor
param.Cg = 0.9;                      % fatigue gradient factor
param.Cs = 0.7;                      % fatigue surface factor
param.primaryHeading = 315;          % primary wind direction [deg]
param.secondaryHeading = 157.5;      % secondary wind direction [deg]

turbine.alpha = 1/3;                      % axial induction factor
turbine.nBlades = 3;                      % number of turbine blades

turbine.maxRPM = 15.5;                    % max turbine speed [RPM]
turbine.minRPM = 9.6;                     % min turbine [RPM]
turbine.P = 2.5;                          % Rated power [MW]

turbine.R = 48 ;                          % Blade Radius [m]
turbine.hubR = 1.3;                       % Hub Radius [m]
turbine.SweptArea = pi * turbine.R^2;     % Rotor swept area [m^2]
turbine.hubH = 80.4;                      % Hub height [m]
turbine.towerH = 77.7;                    % Height of tower portion [m]
turbine.nacelleH = 2*(turbine.hubH - ...  % Height of nacelle [m]
    turbine.towerH);




%% DELIVERABLE 1:

%Number of discrete points for blade analysis
nPoints = 300;

%Generate vector of radial distances r along blade
r = linspace(turbine.hubR, turbine.R, nPoints);

%Add values for first deliverable to turbine struct
turbine.pitch = 0;         % pitch, [deg] for first deliverable
turbine.U = 10;            % wind speed [m/s] for first deliverable
turbine.W = 14*2*pi/60;    % Rotational Speed [Rad/s] for first deliverable

%Calculate Cp and Ct
[Q1.C_p, Q1.C_t] = calcCpCt(r, turbine, param);

%Print Cp and Ct to the Command Window
fprintf("Deliverable 1: \nThe Coefficient of Power is: %.4f \n"  + ...
    "The Coefficient of Thrust is: %.4f \n \n", Q1.C_p, Q1.C_t);


% --- DELIVERABLE 1 PLOTS ------

% PLOT DATA GENERATION

%Calculate AoA, Relative wind Speed, relative wind angle. 
[Q1.AoA, Q1.Urel, Q1.windRelAngle] = calcAoA(r, turbine);

%Calculate Section Twist Angle along blade
Q1.twist = calcTwist(r, turbine);

%Calculate Cd, Cl along blade
[Q1.Cd, Q1.Cl] = calcCoef(r, turbine, param);

%Calculate the Loads along blade
[Q1.fTheta, Q1.fZ, Q1.Lift, Q1.Drag] = calcForces(r, turbine, param);

%Generate Plots:
figure(1)
% --- Subplot 1: Aerodynamic Coefficients ---
subplot(3,1,1)
hold on
plot(r, Q1.Cl, 'b', 'LineWidth', 1.4)
plot(r, Q1.Cd, 'g', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Coefficient Value')
title('Aerodynamic Coefficients Along Blade Span')
legend('C_L','C_D','Location','best')
grid on

% --- Subplot 2: Angles ---
subplot(3,1,2)
hold on
plot(r, Q1.windRelAngle, '-b', 'LineWidth', 1.4)
plot(r, Q1.twist, '-r', 'LineWidth', 1.4)
plot(r, Q1.AoA, '-g', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Angle [deg]')
title('Flow Angles Along Blade Span')
legend('\phi (flow angle)','\theta (twist)','\alpha (AoA)','Location','best')
grid on

% --- Subplot 3: Forces ---
subplot(3,1,3)
hold on
plot(r, Q1.Lift, '-b', 'LineWidth', 1.4)
plot(r, Q1.Drag, '-r', 'LineWidth', 1.4)
plot(r, Q1.fZ, '-g', 'LineWidth', 1.4)
plot(r, Q1.fTheta, '-c', 'LineWidth', 1.4)
hold off
xlabel('Distance from Center of Rotation [m]')
ylabel('Load [N/m]')
title('Aerodynamic Loads Along Blade Span')
legend('Lift','Drag','F_z','F_\Theta')
grid on


%% DELIVERABLE 2:

%Team #24 Specific Parameters:
Q2.U = 8.3;                              % deliverable 2 wind speed [m/s]
Q2.tipSpdRatio = 7.27;                   % deliverable 2 tip speed ratio
Q2.W = Q2.tipSpdRatio*Q2.U/turbine.R;    % new rotational speed [rad/s]

%Update Turbine Struct with this new info
turbine.U = Q2.U;
turbine.W = Q2.W;

%Generate pitch sweep array and preallocate Cp array
Q2.pitch = linspace(-10, 10, nPoints); 
Q2.Cp = zeros(1, length(Q2.pitch)); 

%Perform sweep of pitch values
for ii = 1:length(Q2.pitch)
    %Update pitch value
    turbine.pitch = Q2.pitch(ii);

    %Calculate C_p and store in C_p array
    [Q2.Cp(ii), ~] = calcCpCt(r, turbine, param);
end

%Find Max C_p and corresponding pitch angle
Q2.maxCp = max(Q2.Cp);
Q2.maxCpPitch = Q2.pitch(1,(Q2.Cp == Q2.maxCp));

%Print the results to the command window
fprintf("Deliverable 2: \nFor our conditions (%.2f [m/s] windspeed " + ...
    "& %.2f TSR), the max C_p is %.3f, occuring at a pitch of %.3f\n\n", ...
    Q2.U, Q2.tipSpdRatio, Q2.maxCp, Q2.maxCpPitch);

%Plot Cp vs. Pitch Angle
figure(2)
hold on
plot(Q2.pitch, Q2.Cp, LineWidth=2)
plot(Q2.maxCpPitch, Q2.maxCp, 'or', MarkerFaceColor='r')
text(Q2.maxCpPitch-2.1, Q2.maxCp-0.05, sprintf("Max C_p = %.2f",Q2.maxCp) )
xlabel("Pitch Angle")
ylabel("C_p")
title("Deliverable 2: C_p vs. Pitch Angle")

%% DELIVERABLE 3:
%Team 24 Specific Parameters:
Q3.U = 5; %Deliverable 3 wind speed, [m/s]

%Update Turbine Struct with this new info:
turbine.U = Q3.U;

%Generate pitch sweep
Q3.pitch = linspace(-10 ,10, 20); 

%Calculate minimum and maximum Tip Speed Ratio from min/max RPM 
Q3.tpSpdRMin = turbine.minRPM*2*pi / 60 * turbine.R / Q3.U;
Q3.tpSpdRMax = turbine.maxRPM*2*pi / 60 * turbine.R / Q3.U; 

%Generate tip speed ratio sweep
Q3.tpSpdR = linspace(Q3.tpSpdRMin, Q3.tpSpdRMax, 20); 

%preallocate array for 2d Cp data
Q3.CpArray = zeros(length(Q3.pitch), length(Q3.tpSpdR));

%Perform 2d Sweep across pitch and tip speed ratio:
for ii = 1:length(Q3.pitch)
    %Update pitch angle for this iteration
    turbine.pitch = Q3.pitch(ii);

    for jj = 1:length(Q3.tpSpdR)
        %Update rotational speed of turbine for this iteration
        turbine.W = Q3.tpSpdR(jj) * Q3.U / turbine.R;

        %Calculate and store Cp
        [Q3.CpArray(ii,jj), ~] = calcCpCt(r, turbine, param);

    end
end

%Print Maximum C_p to command window
fprintf("Deliverable 3: \nThe maximum C_p for a wind velocity of %.1f "+ ...
    "[m/s] is: %.3f ", Q3.U, max(Q3.CpArray, [], "all"));
fprintf("(See note in report about model breakdown)\n\n")

%Plot Cp vs. Tip Speed Ratio and Pitch Angle
figure(3)
surf(Q3.tpSpdR, Q3.pitch, Q3.CpArray);
title("Deliverable 3: C_p vs. Tip Speed Ratio and Pitch Angle")
xlabel("Tip Ratio [-]")
ylabel("Pitch Angle [deg]")
zlabel("C_p");
colorbar

%% DELIVERABLE 4:
%Team 24 Specific Parameters:
Q4.U = 14;                       % wind speed, m/s for deliverable 4
Q4.W = turbine.maxRPM * 2*pi/60; % max rotational speed (Rad/s) of turbine

%update turbine struct with new data
turbine.U = Q4.U;
turbine.W = Q4.W;

%Find exact pitch to avoid overpowering turbine:
%initial guess, obtained from plotted data
initialGuess = 4;      

%function for use with fzero
zeroFunction = @(x) calcPower(r, ...
    setfield(turbine, 'pitch', x), param) - turbine.P; 

%calculate required pitch with fzero
Q4.reqPitch = fzero(zeroFunction, initialGuess);       


%--- DELIVERABLE 4 PLOTS ---

%set up pitch sweep and preallocate power vector
Q4.pitchsweep = linspace(-10, 10, 100);
Q4.power = zeros([1, length(Q4.pitchsweep)]);

%Sweep through pitch angles and calculate power:
for ii = 1:length(Q4.pitchsweep)
    %update pitch value
    turbine.pitch = Q4.pitchsweep(ii);

    %calculate power
    Q4.power(ii) = calcPower(r, turbine, param);
    
end

%Plot output power vs. pitch angle
figure(4)
hold on
plot(Q4.pitchsweep, Q4.power, LineWidth=2);
yline(turbine.P, '--', Color='k')
plot(Q4.reqPitch, turbine.P,  'or', MarkerFaceColor='r');

title("Deliverable 4: Power Output vs. Pitch Angle")
xlabel("Pitch Angle [deg]");
ylabel("Power Output [MW]");


%% DELIVERABLE 5:
%Tower Structural Analysis: 
Q5.U = Q4.U; %use same wind speed as deliverable 4
Q5.W = Q4.W; %use same max RPM as deliverable 4
Q5.pitch = Q4.reqPitch; %use pitch found in deliverable 4

%update struct
turbine.U = Q5.U;
turbine.W = Q5.W;
turbine.pitch = Q5.pitch;

%generate h array over height of entire tower and nacelle
h = linspace(1e-6, (turbine.towerH + turbine.nacelleH), nPoints);

%calculate load, shear, moment, angle, and deflection for all h
[Q5.q,Q5.v,Q5.m,Q5.theta,Q5.delta] = calcBendingAnalysis(h,r,turbine,param);

%find max deflection in tower portion only
Q5.maxDeflection = max(abs(Q5.delta(h<=turbine.towerH)));

%calculate stress in tower portion
Q5.towerStress = calcTowerStress(h,r, turbine, param);

%calculate safety factor against fatigue failure
Q5.fatigueSF = calcFatigueSF(h,r,turbine,param);


% % % % I = calcI(h, turbine);
% % % % t = calcTowerT(h(h<=turbine.towerH), turbine);
% % % % figure(6)
% % % % plot(h(h<=turbine.towerH), I(h<=turbine.towerH))
% % % % figure(7)
% % % % hold on
% % % % plot(h(h<=turbine.towerH), t)
% % % % plot(turbine.towerSpecs.("Height (mm)")/1000, turbine.towerSpecs.("Wall thk (mm)")/1000);


% --- DELIVERABLE 5 PLOTS ---
hTower = h(h<=turbine.towerH);

figure(5)
subplot(5,1,1)
plot(hTower, Q5.q(h<=turbine.towerH))
xlabel("Height [m]")
ylabel("q [N/m]")

subplot(5,1,2)
plot(hTower, Q5.v(h<=turbine.towerH));
xlabel("Height [m]")
ylabel("Shear Force, [N]");

subplot(5,1,3)
plot(hTower, Q5.m(h<=turbine.towerH));
xlabel("Height [m]")
ylabel("Bending Moment, [N*m]");

subplot(5,1,4)
plot(hTower, Q5.theta(h<=turbine.towerH));
xlabel("Height [m]")
ylabel("Deflection Angle, [rad]")

subplot(5,1,5)
plot(hTower, Q5.delta(h<=turbine.towerH));
xlabel("Height [m]")
ylabel("Deflection [m]")




%% -------------FUNCTIONS:------------------ %%

%% calcTwist
%calculate airfoil twist for every r along blade

function [sectionTwist] = calcTwist(r, turbine)
    r = r*1000; %convert r to mm
    BladeData = turbine.profileData;
    sectionTwist = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.BladeTwist, r);
end

%% calcChord
%calculate chord length for every r along blade

function [chordLength] = calcChord(r, turbine)
    r = r*1000; %convert r to mm
    BladeData = turbine.profileData;
    chordLength = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.ChordLength, r);
    chordLength = chordLength / 1000; %return chord in meters
end

%% calcAPrime
%calculate angular induction factor along blade
function alphaPrime = calcAPrime(r, turbine)

    W = turbine.W; %turbine rotational speed, rad/s
    lambdaR = (W .*r)/turbine.U; %get local speed ratio
    alpha = turbine.alpha; %get axial induction factor

    alphaPrime = -1/2 + 1/2.*sqrt(1+4./((lambdaR).^2).*alpha.*(1-alpha)); %calculate angular induction factor
  
end

%% CalcAoA
%calculate the AoA for every r along blade

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

%% getAirfoil
%return airfoil type for the given r value(s)

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

%% calcCoef
%calculate the drag and lift coefficients for every r along blade

function [Cd, Cl] = calcCoef(r, turbine, param)
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
            Re = param.rho_air .* Urel(ii) .* chords(ii) / param.mu_air;
            
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

%% calcForces
%calculate the distributed forces [N/m] for every r along blade

function [fTheta, fZ, Lift, Drag] = calcForces(r, turbine, param)
    [~,Urel,windRelAngle]  = calcAoA(r,turbine); % get relative wind speed and angle (deg)
    chordLength = calcChord(r,turbine); %get chord length for every r
    [Cd, Cl] = calcCoef(r, turbine, param); %get coefficients
    
    
    Lift = 0.5 .* param.rho_air .* Urel.^2 .* Cl .* chordLength;
    Drag = 0.5 * param.rho_air .* Urel.^2 .* Cd .* chordLength;

    fTheta = Lift .* sind(windRelAngle) - Drag .* cosd(windRelAngle);
    fZ = Lift .*cosd(windRelAngle) + Drag.* sind(windRelAngle);

end

%% calcPower
%calculate power output of turbine

function power = calcPower(r, turbine, param)
    [fTheta, ~, ~, ~] = calcForces(r, turbine, param);
    
    d_Torque = fTheta .* r;
    Torque = trapz(r, d_Torque);
    power = turbine.nBlades* Torque * turbine.W ./ 10^6; %return power in MW


end

%% calcThrust
%calculate thrust load of turbine

function thrust = calcThrust(r, turbine, param)
    [~, fZ, ~, ~] = calcForces(r, turbine, param);
    thrust = turbine.nBlades* trapz(r, fZ);
end

%% calcCpCt
%calculate coefficient of power and thrust

function [C_p, C_t] = calcCpCt(r, turbine, param)

    %calculate power, convert from [MW] to [W]
    Power = calcPower(r, turbine, param) * 10^6; 

    %calculate thrust
    Thrust = calcThrust(r, turbine, param);

    %calculate coefficients
    C_p = Power / (1/2 * param.rho_air * turbine.U^3 ...
        * turbine.SweptArea);
    C_t = Thrust / (1/2 * param.rho_air * turbine.U^2 ...
        * turbine.SweptArea);
end

%% calcWind
%return the speed of the wind [m/s] at height h [m]

function windSpeed = calcWind(h, turbine)
    refSpeed = turbine.U; %reference wind speed, m/s
    refHeight = turbine.hubH; %reference height (hub height)

    windSpeed = refSpeed .* (h./refHeight).^(1/7);
end

%% calcDia
%return the diameter of the tower in meters for h given in meters

function towerDia = calcTowerDia(hTower, turbine)
    heights = hTower*1000; %convert h vector to mm
    heightData = turbine.towerSpecs.("Height (mm)");
    diaData = turbine.towerSpecs.("OD (mm)");
    towerDia = interp1(heightData, diaData, heights)./1000; %return dia in meters
end

%% calcTowerT
%calculate thickness for h along tower portion

function towerT = calcTowerT(hTower, turbine)
    heights = hTower*1000;
    heightData = turbine.towerSpecs.("Height (mm)");
    tData = turbine.towerSpecs.("Wall thk (mm)");
    towerT = interp1(heightData, tData, heights)./1000; %return tower thickness in meters
end

%% CalcTowerCd
%calculate drag coefficient for h along tower portion

function Cd = calcTowerCd(hTower, turbine, param)
    windSpeeds = calcWind(hTower, turbine); %calculate wind speed in m/s

    dia = calcTowerDia(hTower, turbine); %calculate diameters along tower h

    Re = param.rho_air .* windSpeeds .*dia / param.mu_air; %calculate Re along tower h

    Cd = cylinderCD(Re); %calculate Cd
end

%% calcTowerLoad
%returns the load on the tower for all h, including nacelle portion

function towerLoad = calcTowerLoad(h, r, turbine, param)

    %get Cd, diameter, wind speed along tower portion
    hTower = h(h<=turbine.towerH);
    Cd = calcTowerCd(hTower, turbine, param ...
        );
    towerDia = calcTowerDia(hTower, turbine);
    windSpeeds = calcWind(hTower, turbine);
    
    %calculate drag load [n/m] along tower portion
    towerDrag = 0.5 * param.rho_air .* windSpeeds.^2 .* Cd .* towerDia;

    %calculate the thrust load of the nacelle
    totalThrust = calcThrust(r, turbine, param);
    distThrust = totalThrust / turbine.nacelleH;
    
    %add distributed nacelle load to drag load array
    towerLoad = paddata(towerDrag, length(h), FillValue=distThrust);

end

%% CalcI
%returns the bending moment [m^4] for all h, including nacelle portion
function I = calcI(h, turbine)

    %calculate diameter and thickness of tower portion (excl. nacelle)
    hTower = h(h<=turbine.towerH);
    diameters = calcTowerDia(hTower, turbine);
    thickness = calcTowerT(hTower, turbine);
    
    %calculate bending moment for tower portion
    TowerI = pi.*(diameters/2).^3 .* thickness;

    %define large I value for nacelle
    nacelleI = 1e25; 
    
    %append nacelle values to end of I array
    I = paddata(TowerI, length(h), FillValue=nacelleI);

end

%-----------
%SHEAR, MOMENT, THETA, DEFLECTION FUNCTIONS HERE (DONT FORGET TO MULTIPLY I
%by E)

%% calcBendingAnalysis
%calculate the shear, moment, angle, and deflection along the tower

function [q,v,m,theta,delta] = calcBendingAnalysis(h, r, turbine, param)

    %calculate distributed load on tower 
    q = -calcTowerLoad(h, r, turbine, param);

    %calculate flexural stiffness
    EI = calcI(h, turbine).*param.E_steel*10^9;

    %Find Shear as a function of h, dV/dx=-q -> V=-q*x+C1
    %reaction (initial) shear is the sum of all acting forces (Drag+Thrust)
    initialV=-trapz(h,q);
    
    %Now shear as a function of h, each value is the internal shear at H
    v=cumtrapz(h,q)+initialV;
    
    %Find Moment as a function h, dM/dx=V -> M=V*x+C1*x+C2
    %initial condition of moment is sum of all individual moments
    initialM=-trapz(h,q.*h);
    
    %calculate M by cummulative integration of v along tower
    m = cumtrapz(h,v)-initialM;

    %calculate deflection angle along h
    theta = cumtrapz(h, m) ./EI;

    %calculate deflection along h
    delta = cumtrapz(h, theta);

end

%-----------

%% calcTowerStress
%calculate the maximum stress seen for every h along the tower
function towerStress = calcTowerStress(h,r,turbine, param)

    %only use portion of h along tower for calculations
    hTower = h(h<=turbine.towerH);

    %calculate moments as function of h (will need h, r, turbine,
    %parameters)
    [~,m,~,~] = calcBendingAnalysis(h,r,turbine,param);

    %take only the moments along the tower
    mTower = m((h<=turbine.towerH));

    %calculate the distance from neutral axis to outer surface, c
    c = calcTowerDia(hTower, turbine)/2; 

    I = calcI(hTower, turbine);

    towerStress = mTower .* c ./ I;

end

%% calcFatigueS
%calculate mean and alternating stresses on the tower

function [meanS, altS] = calcFatigueS(h,r,turbine,param)
    %calculate max (tensile) stress 
    maxStress = max(calcTowerStress(h, r, turbine, param));

    %calculate min stress from angular difference in wind headings
    minStress = maxStress*cosd(param.primaryHeading-param.secondaryHeading);

    %calculate the mean stress
    meanS = maxStress + minStress / 2;

    %calculate the alternating stress
    altS = abs(maxStress - meanS);

end

%% calcFatigueSF
%calculate the safety factor against fatigue failure

function fatigueSF = calcFatigueSF(h, r, turbine, param)
    %calculate mean and alternating stress
    [meanS, altS] = calcFatigueS(h, r, turbine, param);

    %calculate endurance strength of material
    sigma_e = param.TS_steel*param.Cg*param.Cl*param.Cs*10^6;

    %retrieve tensile strength of material
    TS = param.TS_steel*10^6;

    %calculate Safety Factor
    fatigueSF = 1/( (altS/sigma_e) + (meanS/TS) );

end