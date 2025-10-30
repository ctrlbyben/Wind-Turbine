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
param.steelYield = 345;            % yield strength of steel [MPa]
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
nPoints = 101;

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


% % % % % --- DELIVERABLE 1 PLOTS ------
% % % % 
% % % % % PLOT DATA GENERATION
% % % % 
% % % % %Calculate AoA, Relative wind Speed, relative wind angle. 
% % % % [Q1.AoA, Q1.Urel, Q1.windRelAngle] = calcAoA(r, turbine);
% % % % 
% % % % %Calculate Section Twist Angle along blade
% % % % Q1.twist = calcTwist(r, turbine);
% % % % 
% % % % %Calculate Cd, Cl along blade
% % % % [Q1.Cd, Q1.Cl] = calcCoef(r, turbine, param);
% % % % 
% % % % %Calculate the Loads along blade
% % % % [Q1.fTheta, Q1.fZ, Q1.Lift, Q1.Drag] = calcForces(r, turbine, param);
% % % % 
% % % % %Generate Plots:
% % % % figure(1)
% % % % % --- Subplot 1: Aerodynamic Coefficients ---
% % % % subplot(3,1,1)
% % % % hold on
% % % % plot(r, Q1.Cl, 'b', 'LineWidth', 1.4)
% % % % plot(r, Q1.Cd, 'g', 'LineWidth', 1.4)
% % % % hold off
% % % % xlabel('Distance from Center of Rotation [m]')
% % % % ylabel('Coefficient Value')
% % % % title('Aerodynamic Coefficients Along Blade Span')
% % % % legend('C_L','C_D','Location','best')
% % % % grid on
% % % % 
% % % % % --- Subplot 2: Angles ---
% % % % subplot(3,1,2)
% % % % hold on
% % % % plot(r, Q1.windRelAngle, '-b', 'LineWidth', 1.4)
% % % % plot(r, Q1.twist, '-r', 'LineWidth', 1.4)
% % % % plot(r, Q1.AoA, '-g', 'LineWidth', 1.4)
% % % % hold off
% % % % xlabel('Distance from Center of Rotation [m]')
% % % % ylabel('Angle [deg]')
% % % % title('Flow Angles Along Blade Span')
% % % % legend('\phi (flow angle)','\theta (twist)','\alpha (AoA)','Location','best')
% % % % grid on
% % % % 
% % % % % --- Subplot 3: Forces ---
% % % % subplot(3,1,3)
% % % % hold on
% % % % plot(r, Q1.Lift, '-b', 'LineWidth', 1.4)
% % % % plot(r, Q1.Drag, '-r', 'LineWidth', 1.4)
% % % % plot(r, Q1.fZ, '-g', 'LineWidth', 1.4)
% % % % plot(r, Q1.fTheta, '-c', 'LineWidth', 1.4)
% % % % hold off
% % % % xlabel('Distance from Center of Rotation [m]')
% % % % ylabel('Load [N/m]')
% % % % title('Aerodynamic Loads Along Blade Span')
% % % % legend('Lift','Drag','F_z','F_\Theta')
% % % % grid on


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
    "& %.2f TSR), the max C_p is %.3f, occuring at a pitch of %.3f " + ...
    "[deg] \n\n", Q2.U, Q2.tipSpdRatio, Q2.maxCp, Q2.maxCpPitch);

%--- DELIVERABLE 2 PLOTS

%Plot Cp vs Pitch
figure(2)
hold on
plot(Q2.pitch, Q2.Cp, LineWidth=2)
plot(Q2.maxCpPitch, Q2.maxCp, 'or', MarkerFaceColor='r', MarkerSize=8)
text(Q2.maxCpPitch-2.1, Q2.maxCp-0.05, sprintf("Max C_p = %.2f", ...
    Q2.maxCp), FontName='Times New Roman', FontSize=10)

%plot formatting and auto-save
xlabel("Pitch Angle [deg]")
ylabel("C_p")
legend(["Coefficient of Power (C_p", "Maximum Value"], Location= ...
    'southwest')
set(gca,'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 4];
print(fig, 'Cp_Pitch.png', '-dpng', '-r300');


%% DELIVERABLE 3:
%Team 24 Specific Parameters for deliverable 4:
Q3.U = 5; %wind speed, [m/s]

%Update Turbine Struct:
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
fprintf("Deliverable 3: \nThe maximum C_p for a wind velocity of %.1f "+...
    "[m/s] is: %.3f ", Q3.U, max(Q3.CpArray, [], "all"));
fprintf("(See note in report about model breakdown)\n\n")

% --- DELIVERABLE 3 PLOTS---

%Plot Cp vs. Tip Speed Ratio and Pitch Angle
Q3.CpArray(Q3.CpArray<0) = NaN; %remove negative values
figure(3)
surf(Q3.tpSpdR, Q3.pitch, Q3.CpArray);
view([-37.5-100, 30])

%plot formatting and auto-save
xlabel("Tip Ratio [-]")
ylabel("Pitch Angle [deg]")
zlabel("C_p");
colorbar
set(gca,'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 4];
print(fig, 'Cp_TSR_pitch', '-dpng', '-r300');

%% DELIVERABLE 4:
%Team 24 Specific Parameters:
Q4.U = 14;                       % wind speed, m/s for deliverable 4
Q4.W = turbine.maxRPM * 2*pi/60; % max rotational speed (Rad/s) of turbine

%update turbine struct with new data
turbine.U = Q4.U;
turbine.W = Q4.W;

%Find exact pitch to avoid overpowering turbine:
%initial guess, obtained from plotted data
initialGuess1 = 4;  
initialGuess2 = -6;

%function for use with fzero
zeroFunction = @(x) calcPower(r, ...
    setfield(turbine, 'pitch', x), param) - turbine.P; 

%calculate required pitch with fzero
Q4.reqPitch1 = fzero(zeroFunction, initialGuess1);     
Q4.reqPitch2 = fzero(zeroFunction, initialGuess2);

%Print Result to command window
fprintf("Deliverable 4:\nFor a windspeed of %.1f [m/s] at max RPM, " + ...
    "the blades must be pitched to %.3f [deg] or %.3f [deg] to avoid " +...
    "overpowering the turbine. \n\n",Q4.U, Q4.reqPitch1, Q4.reqPitch2);

%--- DELIVERABLE 4 PLOTS ---

%GENERATE PITCH SWEEP DATA FOR PLOTTING:
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
yline(turbine.P, '--', Color='k', LineWidth =2, Layer='bottom')
plot([Q4.reqPitch1,Q4.reqPitch2] , [turbine.P, turbine.P], ...
    'or', MarkerFaceColor='r', MarkerSize=8);

%formatting and auto-save
xlabel("Pitch Angle [deg]");
ylabel("Power Output [MW]");
legend(["Turbine Power", "Rated Power", "Pitch Angle Limits"], ...
    Location="southwest")
set(gca,'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 4];
print(fig, 'pwr_pitch', '-dpng', '-r300');

%% DELIVERABLE 5:
%Tower Structural Analysis: 
Q5.U = Q4.U; %use same wind speed as deliverable 4
Q5.W = Q4.W; %use same max RPM as deliverable 4
Q5.pitch = Q4.reqPitch1; %use positive pitch found in deliverable 4

%update struct with new values for deliverable 5:
turbine.U = Q5.U;
turbine.W = Q5.W;
turbine.pitch = Q5.pitch;

%generate h array over height of entire tower and nacelle
h = linspace(1e-6, (turbine.towerH + turbine.nacelleH), nPoints);

%calculate load, shear, moment, angle, and deflection for all h
[Q5.q,Q5.v,Q5.m,Q5.theta,Q5.delta] = calcBendingAnalysis(h,r,turbine, ...
    param);

%find max deflection in tower portion only
Q5.maxDeflection = max(abs(Q5.delta(h<=turbine.towerH)));

%calculate stress in tower portion
Q5.towerStress = calcTowerStress(h,r, turbine, param);

%calculate max stress in tower
Q5.maxStress = max(Q5.towerStress);

%calculate safety factor against yielding
Q5.SF = param.steelYield * 10^6 / Q5.maxStress;

%calculate safety factor against fatigue failure
Q5.fatigueSF = calcFatigueSF(h,r,turbine,param);

%Print results to command window:
fprintf("Deliverable 5: \nFor a windspeed of %.1f [m/s] at max RPM," + ...
   " the tower structural performance is: \nMax Deflection: %.3f [m]"+...
    "\nMax Stress: %.3f [MPa] \nSafety Factor against yielding: " + ...
    "%.2f \nSafety Factor in Fatigue: %.2f\n\n", Q5.U, Q5.maxDeflection,...
    Q5.maxStress/10^6, Q5.SF, Q5.fatigueSF);


% --- DELIVERABLE 5 PLOTS ---

%generate data for Goodman Diagram:
[Q5.sigmaMean, Q5.sigmaAlt] = calcFatigueS(h,r,turbine,param);
Q5.sigma_e = 0.5*param.TS_steel*param.Cg*param.Cl*param.Cs;

%Plot Goodman Diagram
figure(5)
hold on
plot([0,param.TS_steel], [Q5.sigma_e, 0], LineWidth=2)
plot(Q5.sigmaMean/10^6, Q5.sigmaAlt/10^6, 'or', MarkerFaceColor='r');
xlabel("\sigma_m [MPa]")
ylabel("\sigma_a [MPa]")
legend(["Line of Constant Life", "Load Point"])

%ormat and auto-save goodman plot
set(gca,'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 4];
print(fig, 'goodman', '-dpng', '-r300');

%Plot load, shear, moment, angle, and deflection over ~only~ tower portion
hTower = h(h<=turbine.towerH); %Tower portion of h array

figure(6)
tiledlayout("vertical")
nexttile
plot(hTower, Q5.q(h<=turbine.towerH), LineWidth=2)
ylabel("q [N/m]")
xlim([0,max(hTower)])

nexttile
plot(hTower, Q5.v(h<=turbine.towerH), LineWidth=2);
ylabel("V [N]");
xlim([0,max(hTower)])

nexttile
plot(hTower, Q5.m(h<=turbine.towerH), LineWidth=2);
ylabel("M [Nm]");
xlim([0,max(hTower)])

nexttile
plot(hTower, Q5.theta(h<=turbine.towerH), LineWidth=2);
ylabel("\theta [rad]")
xlim([0,max(hTower)])

nexttile
plot(hTower, Q5.delta(h<=turbine.towerH), LineWidth=2);
xlabel("Height [m]")
ylabel("\delta [m]")
xlim([0,max(hTower)])

%format bending analysis plots
fig = gcf;
set(findobj(fig,'type','axes'),'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 6];
print(fig, 'bending_analysis', '-dpng', '-r300');

%Plot stress along tower section
figure(7)
hold on
plot(hTower, Q5.towerStress./10^6, LineWidth=2)
plot(hTower(Q5.towerStress == Q5.maxStress), Q5.maxStress/10^6, 'or', ...
    MarkerFaceColor='r', MarkerSize=8);
xlabel("Height [m]")
ylabel("Stress [MPa]")
legend(["Stress in Tower", "Maximum Stress Value"], Location="southwest")
xlim([0,max(hTower)])

%format tower stress plot
set(gca,'FontName','Times New Roman','FontSize',10)
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPositionMode = 'manual';
fig.PaperPosition = [0, 0, 5.5, 2.5];
print(fig, 'stress_analysis', '-dpng', '-r300');


%% -------------FUNCTIONS:------------------ %%

%% calcTwist
function [sectionTwist] = calcTwist(r, turbine)
%CALCTWIST calculates airfoil twist for every r along the blade
%[sectionTwist] = calcTwist(r, turbine) calculates the section twist of the
%turbine blade at the given point(s) r in meters. The output is given in
%degrees. 

    %convert r to mm
    r = r*1000;

    %get blade profile data table
    BladeData = turbine.profileData;

    %interpolate the data table for the twist at r
    sectionTwist = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.BladeTwist, r);
end

%% calcChord
function [chordLength] = calcChord(r, turbine)
%CALCCHORD calculates the chord length at every r along the blade
%[chordLength] = calcChord(r, turbine) interpolates the given blade profile
%data to determine the chord length at the given blade radial position(s)
%r, in meters. The chord length(s) are returned in meters. 

    %convert r to mm
    r = r*1000; 

    %get blade profile data table
    BladeData = turbine.profileData;

    %interpolate the data table for the chord length at r
    chordLength = interp1(BladeData.DistanceFromCenterOfRotation,...
        BladeData.ChordLength, r);

    %return chord in meters
    chordLength = chordLength / 1000; 
end

%% calcAPrime
function [alphaPrime] = calcAPrime(r, turbine)
%ALPHAPRIME calculates the angular induction factor
%[alphaPrime] = calcAPrime(r, turbine) takes the given radial position
%vector r in meters, and calculates the angular induction factor alpha 
% prime for all r. The angular induction factor is unitless.  
    
    %get turbine rotational speed, rad/s
    W = turbine.W; 

    %calculate speed ratio for all r
    lambdaR = (W .*r)/turbine.U; 
    
    %get axial induction factor from turbine struct
    alpha = turbine.alpha; 

    %calculate angular induction factor
    alphaPrime = -1/2 + 1/2.*sqrt(1+4./((lambdaR).^2).*alpha.*(1-alpha)); 
  
end

%% CalcAoA
function [AoA,Urel,windRelAngle]  = calcAoA(r,turbine)
%CALCAOA calculates the angle of attack along the blade
%[AoA,Urel,windRelAngle]  = calcAoA(r,turbine) takes a radial position
%vector r in meters and the turbine struct, and calculates the angle of
%attack [deg] and the relative wind velocity [m/s] and angle [deg]. 

    %get turbine rotational speed [rad/s] and wind speed [m/s]
    W = turbine.W; 
    U = turbine.U; 

    %get axial and angular induction factors
    a = turbine.alpha; 
    a_prime = calcAPrime(r,turbine); 

    %calculate section pitch
    sectionPitch = calcTwist(r, turbine) + turbine.pitch; 

    %calculate the angle of relative wind using the components
    tangentComponent = r.*W.*(1+a_prime);
    axialComponent = (U.*(1-a));
    windRelAngle = atand(axialComponent ./ tangentComponent);

    %calculate the angle of attack
    AoA = windRelAngle-sectionPitch;
    
    %calculate magnitude of relative wind velocity
    Urel = sqrt(axialComponent.^2 + tangentComponent.^2);
end

%% getAirfoil
function [airfoil] = getAirfoil(r)
%GETAIRFOIL calculates the airfoil type across r. 
%[airfoil] = getAirfoil(r) takes the radial position vector r [m] and
% returns the airfoil type as a string.

    %convert r to mm
    r = r*1000; 

    %preallocate array
    airfoil = strings(1,length(r)); 

    %filter r value to get airfoil type
    for ii = 1:length(r) 
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
function [Cd, Cl] = calcCoef(r, turbine, param)
%CALCCOEF calculates the drag and lift coefficients along the blade
%[Cd, Cl] = calcCoef(r, turbine, param) takes the radial position vector r
%[m], the turbine struct, and the parameter struct, and outputs the
%coefficients of drag and lift along the blade.

    %get the airfoil type for all r
    airfoils = getAirfoil(r);

    %get the AoA and relative wind velocity for all r
    [AoA,Urel,~] = calcAoA(r, turbine);

    %get chord lengths
    chords = calcChord(r, turbine);

    %Preallocate arrays to store Cl and Cd
    Cl = zeros(1,length(r)); 
    Cd = zeros(1,length(r)); 

    %Calculate Cd and Cl for all r based on airfoil type and AoA
    for ii = 1:length(r)
        %determine airfoil and AoA for current iteration
        localAirfoil = airfoils(ii);
        localAoA = AoA(ii);
        
        %use provided function if the section is a circle
        if localAirfoil == "circle"
            Cl(ii) = 0;
            %calculate reynolds number of section
            Re = param.rho_air .* Urel(ii) .* chords(ii) / param.mu_air;
            
            %calculate CD using external function
            Cd(ii) = cylinderCD(Re); 

        %otherwise, use airfoil data tables
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
function [fTheta, fZ] = calcForces(r, turbine, param)
%CALCFORCES calculates the forces acting on the blade for all r
%[fTheta, fZ] = calcForces(r, turbine, param) takes the radial
%position vector r [m] and the turbine and parameter structs, and
%calculates the turbine blade load [N/m] in both the axial (fZ) and 
%tangiential (fTheta) direction.

    %get relative wind speed and angle (deg)
    [~,Urel,windRelAngle]  = calcAoA(r,turbine); 

    %get chord length for every r
    chordLength = calcChord(r,turbine);

    %get coefficients
    [Cd, Cl] = calcCoef(r, turbine, param); 
    
    %calculate lift and drag force on the airfoil
    Lift = 0.5 .* param.rho_air .* Urel.^2 .* Cl .* chordLength;
    Drag = 0.5 * param.rho_air .* Urel.^2 .* Cd .* chordLength;

    %calculate the portion of the lift and drag force acting in the
    %tangiential and axial directions. 
    fTheta = Lift .* sind(windRelAngle) - Drag .* cosd(windRelAngle);
    fZ = Lift .*cosd(windRelAngle) + Drag.* sind(windRelAngle);

end

%% calcPower
function [power] = calcPower(r, turbine, param)
%CALCPOWER calculates the power output of the turbine. 
%[power] = calcPower(r, turbine, param) takes the radial position vector r
%[m] and the turbine and parameter structs, and calculates the power output
%of the turbine in MegaWatts. 

    %calculate the load [N/m] in the tangiential direction
    [fTheta, ~] = calcForces(r, turbine, param);
    
    %calculate the differential torque component
    d_Torque = fTheta .* r;

    %integrate differential torque for total torque
    Torque = trapz(r, d_Torque);

    %calculate power with torque, rotational speed, and the number of 
    %blades
    power = turbine.nBlades* Torque * turbine.W ./ 10^6; %[MW]


end

%% calcThrust
function [thrust] = calcThrust(r, turbine, param)
%CALCTHRUST calculates the total thrust force acting on the wind turbine.
%[thrust] = calcThrust(r, turbine, param) takes the radial position vector 
%r [m] along with the turbine and parameter structs, and integrates the 
%distributed axial load to find the total rotor thrust [N].

    %calculate the load [N/m] on the blade in the axial direction
    [~, fZ] = calcForces(r, turbine, param);

    %calculate the total force [N] by integrating across all 3 blades
    thrust = turbine.nBlades* trapz(r, fZ);
end

%% calcCpCt
function [C_p, C_t] = calcCpCt(r, turbine, param)
%CALCCPCT calculates the coefficients of power and thrust for the turbine.
%[C_p, C_t] = calcCpCt(r, turbine, param) uses the turbine geometry, 
%aerodynamic loading, and flow properties to determine the nondimensional 
%power coefficient C_p and thrust coefficient C_t.

    %calculate power, convert from [MW] to [W]
    Power = calcPower(r, turbine, param) * 10^6; 

    %calculate thrust force
    Thrust = calcThrust(r, turbine, param);

    %calculate coefficients
    C_p = Power / (1/2 * param.rho_air * turbine.U^3 ...
        * turbine.SweptArea);
    C_t = Thrust / (1/2 * param.rho_air * turbine.U^2 ...
        * turbine.SweptArea);
end

%% calcWind
function [windSpeed] = calcWind(h, turbine)
%CALCWIND calculates the wind speed at a given height.
%[windSpeed] = calcWind(h, turbine) takes the height h [m] and turbine 
%struct containing the reference wind velocity and hub height, and 
%returns the wind speed [m/s] at that elevation.

    %get reference wind speed and height, [m/s], [m]
    refSpeed = turbine.U; 
    refHeight = turbine.hubH; 

    %calculate wind speed for all h
    windSpeed = refSpeed .* (h./refHeight).^(1/7);
end

%% calcDia
function [towerDia] = calcTowerDia(hTower, turbine)
%CALCTOWERDIA calculates the tower outer diameter as a function of height.
%[towerDia] = calcTowerDia(hTower, turbine) takes the tower height hTower
%[m] and turbine struct containing tower geometry data, and interpolates 
%the outer diameter [m] at each height.

    %Convert h vector to mm
    heights = hTower*1000;

    %get tower section data
    heightData = turbine.towerSpecs.("Height (mm)");
    diaData = turbine.towerSpecs.("OD (mm)");

    %interpolate section data for diameter [m]
    towerDia = interp1(heightData, diaData, heights)./1000; 
end

%% calcTowerT
function towerT = calcTowerT(hTower, turbine)
%CALCTOWERT calculates the tower wall thickness along its height.
%towerT = calcTowerT(hTower, turbine) takes the tower height hTower [m] and 
%the turbine struct containing geometry data, and interpolates the wall 
%thickness [m] for each specified height.

    %convert heights to mm
    heights = hTower*1000;

    %get tower section data
    heightData = turbine.towerSpecs.("Height (mm)");
    tData = turbine.towerSpecs.("Wall thk (mm)");

    %interpolate section data for thickness [m]
    towerT = interp1(heightData, tData, heights)./1000; 
end

%% CalcTowerCd
function Cd = calcTowerCd(hTower, turbine, param)
%CALCTOWERCD determines the drag coefficient of the tower along its height.
%Cd = calcTowerCd(hTower, turbine, param) takes the tower height hTower [m], 
%along with turbine and parameter structs, and computes the local drag 
%coefficient using the provided function

    %calculate wind speed in m/s
    windSpeeds = calcWind(hTower, turbine);

    %calculate diameters along tower h
    dia = calcTowerDia(hTower, turbine); 

    %calculate Re along tower h
    Re = param.rho_air .* windSpeeds .*dia / param.mu_air; 

    %calculate Cd
    Cd = cylinderCD(Re);
end

%% calcTowerLoad
function towerLoad = calcTowerLoad(h, r, turbine, param)
%CALCTOWERLOAD calculates the distributed aerodynamic load along the tower.
%towerLoad = calcTowerLoad(h, r, turbine, param) takes the tower height h 
%[m], radial position vector r [m], and the turbine and parameter structs,
%and computes the drag load [N/m] acting on the tower.

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
function I = calcI(h, turbine)
%CALCI calculates the section bending moment of inertia for all h.
%I = calcI(h, turbine) takes the tower height h[m] and the turbine struct,
%and calculates the bending moment of inertia [m^4] at all positions, 
% including the tower portion.

    %calculate diameter and thickness of tower portion (excl. nacelle)
    hTower = h(h<=turbine.towerH);
    diameters = calcTowerDia(hTower, turbine);
    thickness = calcTowerT(hTower, turbine);
    
    %calculate bending moment for tower portion
    TowerI = pi.*(diameters/2).^3 .* thickness;

    %define very large I value for nacelle
    nacelleI = 1e25; 
    
    %append nacelle values to end of I array
    I = paddata(TowerI, length(h), FillValue=nacelleI); %[m^4]

end

%% calcBendingAnalysis
function [q,v,m,theta,delta] = calcBendingAnalysis(h, r, turbine, param)
%CALCBENDINGANALYSIS performs a deform analysis on the tower
%[q,v,m,theta,delta] = calcBendingAnalysis(h, r, turbine, param) takes in
%the height vector [m], the blade radial vector [r], and the turbine and
%parameter struct, and performs a full deform analysis on the turbine
%tower. The function returns the distributed load q [N/m], the shear force
%v[N], the bending moment m [Nm], the deflection angle theta [rad] and the
%deflection [m]. 


    %calculate distributed load on tower 
    q = -calcTowerLoad(h, r, turbine, param);

    %calculate flexural stiffness
    EI = calcI(h, turbine).*param.E_steel*10^9;

    %reaction (initial) shear is the sum of all acting forces (Drag+Thrust)
    initialV=-trapz(h,q);
    
    %Now shear as a function of h
    v=cumtrapz(h,q)+initialV;
    
    %initial (reaction) moment is sum of all individual moments
    initialM=-trapz(h,q.*h);
    
    %calculate M by cummulative integration of v along tower
    m = cumtrapz(h,v)-initialM;

    %calculate deflection angle along h
    theta = cumtrapz(h, m./EI);

    %calculate deflection along h
    delta = cumtrapz(h, theta);

end

%-----------

%% calcTowerStress
function towerStress = calcTowerStress(h,r,turbine, param)
%TOWERSTRESS calculates the stress along the tower
%towerStress = calcTowerStress(h,r,turbine, param) takes the height and
%radial position vectors, the turbine and parameter structs, and returns
%the maximum stress [Pa] inside the tower cross-section for every height. 

    %only use portion of h along tower for calculations
    hTower = h(h<=turbine.towerH);

    %calculate moments as function of h 
    [~,~,m,~,~] = calcBendingAnalysis(h,r,turbine,param);

    %take only the moments along the tower (excl. nacelle)
    mTower = m((h<=turbine.towerH));

    %calculate the distance from neutral axis to outer surface, c
    c = calcTowerDia(hTower, turbine)/2; 

    %calculate the bending moment of inertia for all h
    I = calcI(hTower, turbine);

    %calculate the stress in the tower for all h
    towerStress = abs(mTower .* c ./ I);

end

%% calcFatigueS
function [meanS, altS] = calcFatigueS(h,r,turbine,param)
%CALCFATIGUES calculates the mean and alternating stresses in the tower
%[meanS, altS] = calcFatigueS(h,r,turbine,param) takes the blade and height
%arrays, and the turbine and parameter structs, and returns the mean and
%alternating stresses in the tower in [Pa]. This calculation involves the
%two primary wind directions, stored inside the turbine struct. 

    %calculate max (tensile) stress 
    maxStress = max(calcTowerStress(h, r, turbine, param));

    %calculate compressive stress from angular difference in wind headings
    minStress = maxStress*cosd(param.primaryHeading- ...
        param.secondaryHeading);

    %calculate the mean stress [Pa]
    meanS = (maxStress + minStress) / 2;

    %calculate the alternating stress [Pa]
    altS = abs(maxStress - meanS)/2;

end

%% calcFatigueSF
function fatigueSF = calcFatigueSF(h, r, turbine, param)
%CALCFATIGUESF calculates the safety factor against fatigue
%fatigueSF = calcFatigueSF(h, r, turbine, param) takes the blade and height
%arrays, and the turbine and parameter struct, and calculates the safety
%factor of the tower against fatigue failure from the changing loads. 

    %calculate mean and alternating stress
    [meanS, altS] = calcFatigueS(h, r, turbine, param);

    %calculate endurance strength of material
    sigma_e = 0.5*param.TS_steel*param.Cg*param.Cl*param.Cs*10^6;

    %retrieve tensile strength of material
    TS = param.TS_steel*10^6;

    %calculate Safety Factor
    fatigueSF = 1/( (altS/sigma_e) + (meanS/TS) );

end