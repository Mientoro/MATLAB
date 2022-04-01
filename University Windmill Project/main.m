clc, clear, close all, format short, format compact

% Peter Erasmus, 30054537
% ENME 337 Final Project
% 

% use xlsread() to read files
% months and airfoil.xlsx

%% Reading Data

% initialize by reading data from files.
disp("Initializing program...");
chord = ReadData("chord.dat", 1);
omega = ReadData("omega.dat", 1);
radius = ReadData("radius.dat", 1);
theta = ReadData("twist.dat", 1);
jan = transpose(xlsread("january.xlsx", "N17:N800"));
feb = transpose(xlsread("february.xlsx", "N17:N800"));
mar = transpose(xlsread("march.xlsx", "N17:N800"));
apr = transpose(xlsread("april.xlsx", "N17:N800"));
may = transpose(xlsread("may.xlsx", "N17:N800"));
june = transpose(xlsread("june.xlsx", "N17:N800"));
july = transpose(xlsread("july.xlsx", "N17:N800"));
aug = transpose(xlsread("august.xlsx", "N17:N800"));
sep = transpose(xlsread("september.xlsx", "N17:N800"));
oct = transpose(xlsread("october.xlsx", "N17:N800"));
nov = transpose(xlsread("november.xlsx", "N17:N800"));
dec = transpose(xlsread("december.xlsx", "N17:N800"));
yearwinds = [jan feb mar apr may june july aug sep oct nov dec];
%airfoil profiles, columns: 1=angle of attack, 2=lift coefficient, 3=drag
du21 = readmatrix("DU21.dat");
du30 = readmatrix("DU30.dat");
du35 = readmatrix("DU35.dat");
du40 = readmatrix("DU40.dat");
naca64 = readmatrix("NACA64.dat");
airfoil = xlsread("Airfoil.xlsx");

disp("Data Retrieved.");

% constants
AIR_DENSITY = 1.23;
CRITICAL_A = 0.2;
CUT_IN_WIND = 3;
CUT_OUT_WIND = 25;
MAX_POWER = 5000000;
VERTICAL_OFFSET = 30;



%% obtain user inputs
hubHeight = input("What will be the hub height in meters? ");
% ensure we have useable input
while hubHeight<=0
    userInput = input("The hub height requires a number. What will be the hub height in meters? ");
    if ~isnan(userInput)
        hubHeight = userInput;
    end
end
blades = input("How many blades will the wind turbine have? ");
% ensure we have useable input
while blades<=0
    userInput = input("The number of blades requires a number. What will be the number of blades? ");
    if ~isnan(userInput)
        blades = userInput;
    end
end
profileNum = menu("Which Airfoil will be used?", "DU21", "DU30", "DU35", "DU40", "NACA64");
profile = du21;
prof = "";
switch profileNum
    case 1
        profile = du21;
        prof = "DU21";
    case 2
        profile = du30;
        prof = "DU30";
    case 3
        profile = du35;
        prof = "DU35";
    case 4
        profile = du40;
        prof = "DU40";
    case 5
        profile = naca64;
        prof = "NACA64";
end


%% Main Code
pAtSpeed = zeros(1,25);

%power of one turbine from windspeeds 1-25 m/s from users selection
for i = 1:length(pAtSpeed)
    pAtSpeed(i) = turbinePower(i, omega(i), radius, chord, theta, blades, profile);
end

%determine yearly power of turbine from user's selection
yearPower = 1:length(yearwinds);
for i = 1:length(yearwinds)
    if isnumeric(yearwinds(i)) && ~isnan(yearwinds(i))
        vHub = yearwinds(i) * (hubHeight/30)^(1/7);
        %convert to m/s
        vHub = vHub / 3.6;
        %ensure windspeed is within limits
        if vHub >= 3 && vHub <= 25
            %add power at this speed
            x = [floor(vHub) ceil(vHub)];
            xq = [pAtSpeed(int32(x(1))), pAtSpeed(int32(x(2)))];
            yearPower(i) = interp1(x, xq, vHub);
        else
            yearPower(i) = 0;
        end
    else
        yearPower(i) = 0;
    end
end

%determine Watt Hours for a single turbine
yearHours = 365 * 24;
%determine # of turbines required to produce power for city for a year with
%city using 16.5 MWh
turbinesNeeded = ceil(16500000*yearHours/sum(yearPower));


%% Visualizations

%output user selection and number of turbines required
fprintf("\nThe selected airfoil is %s.\n", prof);
fprintf("The hub height of the turbine is %.2f m.\n", hubHeight);
fprintf("The turbine has %i blades.\n", blades); 
fprintf("The number of turbines required to sustain a city's power is %i.\n\n", turbinesNeeded);


subplot(2,2,1);
plot(profile(:,1), profile(:,2), profile(:,1), profile(:,3));
title("Lift/Drag Coefficients vs Angle of Attack");
xlabel("Angle of Attack (Degrees)");
ylabel("Lift/Drag Coefficient");
legend("Lift", "Drag");
axis([-200 200 -2.5 2.5]);


subplot(2,2,2);
plot(1:25, pAtSpeed/1000000);
title("Power of Turbine at Wind Speed");
xlabel("Wind Speed (m/s)");
ylabel("Power Produced (MW)");
axis([0 26 0 6]);


subplot(2,2,3);
yyaxis left
histogram(yearwinds, 70);
title("Relative Distribution of Wind Speeds");
xlabel("Wind Speed (km/h)");
ylabel("Occurences");
yyaxis right
scatter(yearwinds, yearPower/1000000, 1, 'filled');
ylabel("Energy (MJ)");
ylim([0 6]);


subplot(2,2,4);
plot(airfoil(:,(2*profileNum-1)), airfoil(:,(2*profileNum)));
title("Airfoil Profile");
xlabel("X");
ylabel("Y");
axis([-0.2 1.2 -0.25 0.25]);





%% Functions


function [power] = turbinePower(v, w, r, c, ang, b, profile)
    pMAX = 5000000; %max generated power
    power = 0;
    if v >= 3 && v <= 25
        w = w * 2 * pi / 60; %convert to radians/second
        forces = zeros(1, length(r));
        torque = zeros(1, length(r)-1);
        %calculate tangential forces
        for i=1:length(r)
            forces(i) = dF(v, w, r(i), c(i), ang(i), b, profile);
        end
        %calculate torque for each segment of the blade
        for i=1:(length(r)-1)
            Ai = (forces(i+1)-forces(i))/(r(i+1)-r(i));
            Bi = (r(i+1)*forces(i)-r(i)*forces(i+1))/(r(i+1)-r(i));
            torque(i) = 1/3*Ai*(r(i+1)^3-r(i)^3)+1/2*Bi*(r(i+1)^2-r(i)^2);
        end
        totalTorque = sum(torque)*b;
        power = totalTorque*w;
        if power > pMAX
            power = pMAX;
        end
    end
end

function [load] = dF(v, w, r, c, ang, b, profile)
    pAIR = 1.23; %air density
    aC = 0.2; %critical a value
    load = 0;
    if v >= 3 && v <= 25
        a = 0.01;
        a1 = 0.01;
        looping = true;
        while looping
            fAng = atand(((1-a)*v)/((1+a1)*w*r));
            attAng = fAng - ang;
            cL = interp1(profile(:,1), profile(:,2), attAng);
            cD = interp1(profile(:,1), profile(:,3), attAng);
            cN = cL*cosd(fAng)+cD*sind(fAng);
            cT = cL*sind(fAng)-cD*cosd(fAng);
            sig = c*b/(2*pi*r);
            A = 1 / (4*sind(fAng)^2/(sig*cN)+1);
            A1 = 1 / (4*sind(fAng)*cosd(fAng)/(sig*cT)-1);
            %correction for large a values
            if (A > aC)
                k = (4*sind(fAng)^2)/(sig*cN);
                A = 1/2 * (2 + k*(1-2*aC) - realsqrt( (k*(1-2*aC)+2)^2 + 4*(k*aC^2-1)));
            end
            % check for values converging
            if (abs(a-A) <= 10^-6) || (abs(a1-A1) <= 10^-6)
                looping = false;
            end
            a = A;
            a1 = A1;
        end
        relV = (v*(1-a))/sind(fAng);
        load = 0.5*cT*pAIR*relV^2*c;
    end
end

