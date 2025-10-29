%% ME5053 PROJECT 2
%Wind Turbine Analysis
%Ben Lahyani, Simon Granberg, Austin Gage

clear; clc; close all;

%% LOAD DATA

turbine.profData = readtable("BladeProfile.csv");
turbine.DU91_W2_250 = readtable("DU91-W2-250.csv");
turbine.DU93_W_210 = readtable("DU93-W-210.csv");
turbine.DU96_W_180 = readtable("DU96-W-180.csv");
turbine.DU97_W_300 = readtable("DU97-W-300.csv");
turbine.towerSpecs = readtable("towerSpecs.csv", VariableNamingRule="preserve"); 

%% CONSTANTS
parameters.rho_air = 1.1;         % Air density [kg/m^3]
parameters.g = 9.81;                % Gravity [m/s^2]
parameters.mu_air = 1.8 * 10^-5;     % [Ns/m^2]
parameters.rho_steel = 7850;      %steel density kg/m^3
parameters.TS_steel = 450; %tensile strength of steel, MPa
parameters.sigmaY_steel = 345;   %yield strength of steel, MPa
parameters.E_steel = 200; %youngs modulus of steel, GPa

turbine.alpha = 1/3;                %axial induction factor
turbine.nBlades = 3;
turbine.R = 48 ; %radius in M
turbine.SweptArea = pi * turbine.R^2;            % Rotor swept area [m^2]
turbine.maxRPM = 15.5; %max turbine speed, RPM
turbine.minRPM = 9.6; %min turbine RPM
turbine.P = 2.5e6;               % Rated power [W]
turbine.H = 80.4;                % Hub height [m]




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
[Q1.AoA, Q1.Urel, Q1.windRelAngle] = calcAoA(r, turbine);

%Calculate Section Twist Angle --------------------------------------------ONLY FOR PLOT Remove for Final 
Q1.twist = calcTwist(r, turbine);

%Calculate Cd, Cl ---------------------------------------------------------ONLY FOR PLOT Remove for Final 
[Q1.Cd, Q1.Cl] = calcCoef(r, turbine, parameters);

%Calculate the Loads on one blade -----------------------------------------ONLY FOR PLOT Remove for Final 
[Q1.fTheta, Q1.fZ, Q1.Lift, Q1.Drag] = calcForces(r, turbine, parameters);

%calculate the power [W] and Thrust Load [N] ------------------------------ONLY FOR PLOT Remove for Final 
[Q1.Power, Q1.Thrust] = calcPower(r, turbine, parameters);

%Calculate Cp and Ct
[Q1.C_p, Q1.C_t] = calcCpCt(r, turbine, parameters);

fprintf("Deliverable 1: \nThe Coefficient of Power is: %.4f \n"  + ...
    "The Coefficient of Thrust is: %.4f \n \n", Q1.C_p, Q1.C_t);

% ----------------Sanity Check Plots-------------------------------------- REMOVE FOR FINAL REPORT
figure(4)
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
    [Q2.C_p, ~] = calcCpCt(r, turbine, parameters);
    Q2.CpSweep(ii) = Q2.C_p;
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
        Q3.currentTipSpeedRatio = Q3.tpRatioSweep(jj);

        %calculate rotational speed, rad/s
        Q3.currentW = Q3.currentTipSpeedRatio * Q3.U / turbine.R;

        %set rotational speed in turbine struct for this iteration
        turbine.W = Q3.currentW;
        [Q3.Cp, ~] = calcCpCt(r, turbine, parameters);

        %store value of Cp in array
        Q3.CpArray(ii,jj) = Q3.Cp;
    end
end

figure(2)
surf(Q3.tpRatioSweep, Q3.pitchSweep, Q3.CpArray);
xlabel("Tip Ratio [-]")
ylabel("Pitch Angle [deg]")
zlabel("C_p");

%% DELIVERABLE 4
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


%% DELIVERABLE 5
%Tower Analysis: 
% Steps: create vector for height positions (h). Then determine wind speed
% at every h. Determine reynolds # for every h, then calculate Cd for each
% h. Determine tower diameter for each h, then calculate drag load (N/m).
% Calculate total thrust load of wind turbine, and assign that to a
% distributed load at the top of the beam. Then, calculate the EI for every
% height (use near infinity for nacelle). 
Q5.h_interval = 0.5;
Q5.U = Q4.U; %use same wind speed as deliverable 4
Q5.W = Q4.W; %use same max RPM as deliverable 4
Q5.pitch = -0.5; %-------------------------------------------------------------------%REPLACE ME ONCE PITCH ANGLE FOR DELIVERABLE 4 IS FOUND

%update struct
turbine.U = Q5.U;
turbine.W = Q5.W;
turbine.pitch = Q5.pitch;

nacelleBottom = (max(turbine.towerSpecs{:,"Height (mm)"})/1000); %determine height of bottom of nacelle
nacelleHeight = 2* (turbine.H - nacelleBottom); %calculate height of nacelle portion
h = linspace(0.000001, (nacelleBottom+nacelleHeight), nPoints);


Q5.windSpeeds = calcWind(h(h<=nacelleBottom), turbine); %calculate wind speed in m/s

Q5.dia = calcTowerDia(h(h<=nacelleBottom), turbine); %calculate diameters along tower h

Q5.Re = parameters.rho_air .* Q5.windSpeeds .* Q5.dia / parameters.mu_air; %calculate Re along tower h

Q5.Cd = calcTowerCd(h(h<=nacelleBottom), turbine, parameters); %determine drag coef along tower portion

Q5.TowerDrag = calcTowerDrag(h(h<=nacelleBottom), turbine, parameters); %determine drag force along tower portion

[Q5.EI, Q5.I]=calcEI(h, turbine, parameters); % determine EI and I values for deflection and bending stress analysis



%add thrust load as distributed load in nacelle portion of h
arrayLength = length(h)-length(Q5.TowerDrag);
[~,Q5.thrustLoad] = calcPower(r, turbine, parameters);
thrustMag = Q5.thrustLoad / nacelleHeight; %determine magnitude of distributed load on nacelle
thrustArray = ones(1,arrayLength) .* thrustMag; %create array for distributed load
TotalTowerDrag = [Q5.TowerDrag, thrustArray];

%calculate EI vector for every h
% NEXT STEPS: use function that calculates EI when in the tower portion,
% and hardcodes a very high value for EI when in the nacelle portion. Will
% need a function to determine wall thickness as well for the EI function. 
%EI DONE (I think, it ends a little short but at the same place as the
%others

%Find Shear as a function of h, dV/dx=-q -> V=-q*x+C1
%Initial condition of shear is the sum of all acting forces (Drag+Thrust)
Q5.initialV=trapz(TotalTowerDrag);

%Now shear as a function of h, each value is the internal shear at H
Q5.V=cumtrapz(-TotalTowerDrag)+Q5.initialV;

%Find Moment as a function h, dM/dx=V -> M=V*x+C1*x+C2
%initial condition of moment is sum of all individual moments
Q5.initialM=trapz(Q5.V.*h);

%Now M as a function of h (cant just multiply v by x, not the same as
%                                                       integration)
Q5.M= cumtrapz(0.5*-TotalTowerDrag.*h.^2+Q5.initialV.*h+Q5.initialM);

%Now can calculate bending stress at each point sigma=M*0.5*OD/I

Q5.sigmaBend = (Q5.M.*calcTowerDia(h, turbine)*0.5)./Q5.I;
Q5.sigmaMax=max(Q5.sigmaBend);


%%Fatigue analysis at our windspeed 
% sigmaMax is sigmamax in primary direction (northwest) lec example 315dg
% sigmaMin is sigmamax in sec direction (use south or SE) lec uses 150 dg
% Will use 165 angle (315-150)
Q5.sigmaAngle= deg2rad(165);

Q5.sigmaMin= Q5.sigmaMax*cos(Q5.sigmaAngle);

Q5.sigmaMean= (Q5.sigmaMax+Q5.sigmaMin)/2;
Q5.sigmaAlt= (Q5.sigmaMax-Q5.sigmaMin)/2;

%Construct Goodman Diagram
%Correction Factors (From lecture)
%Cl=1
%Cg=0.9
%Cs=0.7
Q5.fatigueLimit= 0.5*parameters.TS_steel*10^6*0.9*0.7;

%Make vector to plot goodman line against and plot with load case
x=linspace(0, 450*10^6, 100);

Q5.GoodmanLine=(-Q5.fatigueLimit/(parameters.TS_steel*10^6))*x +Q5.fatigueLimit;

figure(7)
plot(x, Q5.GoodmanLine);
%%UNCOMMENT TO PLOT MARKER
%plot(Q5.sigmaMean, Q5.sigmaAlt, Marker="x", MarkerSize=8)










%--------Sanity Check Plots Deliverable 5-----------
figure(5) 
subplot(5,1,1)
plot(h(h<=nacelleBottom),Q5.windSpeeds)
xlabel("Height [m]")
ylabel("Wind Speed [m/s]")

subplot(5,1,2)
plot(h(h<=nacelleBottom),Q5.dia)
xlabel("Height [m]")
ylabel("Tower Diameter [m]")

subplot(5,1,3)
plot(h(h<=nacelleBottom),Q5.Re)
xlabel("Height [m]")
ylabel("Reynolds Number")

subplot(5,1,4)
plot(h(h<=nacelleBottom), Q5.Cd)
xlabel("Height")
ylabel("C_d")

subplot(5,1,5)
plot(h, TotalTowerDrag)
xlabel("Height")
ylabel("Distributed Load [N/m]")


%EI Sanity Check

    
figure(6)
subplot(4,1,1)
plot(h, Q5.EI)
xlabel("height")
ylabel("EI")

% V sanity check
subplot(4,1,2)
plot(h, Q5.V)
xlabel("height")
ylabel("V(x)")

%M sanity Check
subplot(4,1,3)
plot(h, Q5.M)
xlabel("height")
ylabel("moment")

%Bending Stress sanity check
subplot(4,1,4)
plot(h, Q5.sigmaBend)
xlabel("height")
ylabel("Bending Stress")



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

%% Calculate wind speed for every h along tower section
function windSpeed = calcWind(h, turbine)
    refSpeed = turbine.U; %reference wind speed, m/s
    refHeight = turbine.H; %reference height (hub height)

    windSpeed = refSpeed .* (h./refHeight).^(1/7);
end

%% Calculate Tower Diameter for every h along tower section
function towerDia = calcTowerDia(h, turbine)
    heights = h*1000; %convert h vector to mm
    heightData = turbine.towerSpecs.("Height (mm)");
    diaData = turbine.towerSpecs.("OD (mm)");
    towerDia = interp1(heightData, diaData, heights)./1000; %return dia in meters
end

%% Calculate Tower CD for every h along tower section
function Cd = calcTowerCd(h, turbine, parameters)
    windSpeeds = calcWind(h, turbine); %calculate wind speed in m/s

    dia = calcTowerDia(h, turbine); %calculate diameters along tower h

    Re = parameters.rho_air .* windSpeeds .*dia / parameters.mu_air; %calculate Re along tower h

    Cd = cylinderCD(Re); %calculate Cd
end

%% Calculate Tower Drag Load
function towerDrag = calcTowerDrag(h, turbine, parameters)
    %get Cd, diameter, wind speed for all h
    Cd = calcTowerCd(h, turbine, parameters);
    towerDia = calcTowerDia(h, turbine);
    windSpeeds = calcWind(h, turbine);
    
    %calculate drag load [n/m] on tower for all h
    towerDrag = 0.5 * parameters.rho_air .* windSpeeds.^2 .* Cd .* towerDia;


end

%% Calculate EI for all H
function [EI, I] = calcEI(h, turbine, parameters)

E=parameters.E_steel*10^9;
heights=h*1000;
heightData =turbine.towerSpecs.("Height (mm)");
tdata=turbine.towerSpecs.("Wall thk (mm)");

%Pull Outer Diameters
diamData= calcTowerDia(h, turbine); %in meters
radData=diamData./2; %in meters

thickness= interp1(heightData, tdata, heights)./1000;
innerRads=radData - thickness;%meters

I=(pi/4)*(radData.^4-innerRads.^4); %m^4

EI=E.*I;

end