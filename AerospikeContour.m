% Aerospike Contour Generator

% Constants
M_e = 2.60575;      % Exit Mach Number
gamma = 1.0437;     % Adiabatic Index of the fuel mixture/medium
eta_b = 0;          % Truncation percentile eg. 0.2

% Inputs
A_t = 168;          % Throat Area

% Calculating Area Ratio
AR = Mach2AR(M_e, gamma);

% Calculating Exit Area and Radius
A_e = AR * A_t;
r_e = sqrt(A_e / pi);

% Calculating Prandtl-Meyer angle at the exit
nu_e = Mach2Prandtl(M_e, gamma);

% Number of points for contour generation
N = 25;

% Generating arrays of Mach number, Area Ratio, Prandtl-Meyer angle, and Flow Turn Angle
M_vals = linspace(1, M_e, N);
AR_vals = Mach2AR(M_vals, gamma);
nu_vals = Mach2Prandtl(M_vals, gamma);
mu_vals = Mach2Mangle(M_vals);
alpha_vals = nu_e - nu_vals + mu_vals;

% Calculating non-dimensional values
l_nondim_vals = (1 - sqrt(1 - (AR_vals .* (1 - (eta_b^2)) .* M_vals .* (sin(alpha_vals) ./ AR)))) ./ sin(alpha_vals);
r_nondim_vals = 1 - (l_nondim_vals .* sin(alpha_vals));
x_nondim_vals = l_nondim_vals .* cos(alpha_vals);
y_nondim_vals = l_nondim_vals .* sin(alpha_vals);
Length_nondim = max(x_nondim_vals) - min(x_nondim_vals);

% Calculating dimensional values
l_vals = l_nondim_vals .* r_e;
r_vals = r_nondim_vals .* r_e;
x_vals = x_nondim_vals .* r_e;
y_vals = y_nondim_vals .* r_e;
Length = Length_nondim .* r_e;

% Printing results
fprintf('Exit Mach Number = %g\n', M_e)
fprintf('Area Ratio = %g\n', AR)
fprintf('Length of the spike = %g unit\n', Length)
fprintf('Exit Radius = %g unit\n', r_e)
fprintf('Throat Area = %g sq.units\n', A_t)
fprintf('Exit Area = %g sq.units\n', A_e)
fprintf('Throat Line Length = %g sq.units\n', min(l_vals))
fprintf('Alpha for Throat = %g [deg]\n', alpha_vals(1) * 180 / pi)
fprintf('Flow Turn Angle = %g [deg]\n', nu_e * 180 / pi)

% Create a text file containing coordinates for input in CAD
m = 1;
x_conv = x_vals ./ 10;                   % Value offset for Fusion 360
y_conv = (y_vals - max(y_vals)) ./ 10;   % Value offset for Fusion 360
p = length(x_vals);
x = x_conv(1:m:p);                       % in mm
y = y_conv(1:m:p);                       % in mm
z = zeros(1, N);

% Outputting the x,y,z coordinates of the spike. z is 0. This contour can be revolved in a CAD software to get the 3D profile.
CFDInput = fopen('OG_CFD.csv', 'w');
fprintf(CFDInput, '%6.10f,%12.10f,%12.10f\n', [x; y; z]);
fclose(CFDInput);

%======================= FUNCTIONS =======================%
% Mach number to Prandtl-Meyer angle
function [nu] = Mach2Prandtl(M, gamma) 
    f1 = sqrt((gamma + 1) / (gamma - 1));
    f2 = sqrt((gamma - 1) / (gamma + 1));
    nu = f1 * atan(f2 * sqrt((M.^2) - 1)) - atan(sqrt((M.^2) - 1));
end

% Mach number to Mach angle
function [mu] = Mach2Mangle(M)
    mu = asin(1 ./ M);
end

% Mach number to area ratio (AR)
function [AR] = Mach2AR(M, gamma)
    f1 = 2 / (gamma + 1);
    f2 = (gamma - 1) / 2;
    f3 = (gamma + 1) / (gamma - 1);
    AR = sqrt((1 ./ (M.^2)) .* ((f1 .* (1 + (f2 .* (M.^2)))).^f3));
end

% Area ratio to Mach number
function [M] = AR2Mach(AR, gamma)
    f1 = 2 / (gamma + 1);
    f2 = (gamma - 1) / 2;
    f3 = (gamma + 1) / (gamma - 1);
    M0 = 4;
    M = M0;
    Function = (1 / (M^2)) * ((f1 * (1 + (f2 * (M^2))))^f3) - (AR^2);
    Tolerance = 0.0001;
    maxits = 200;
    J = abs(Function);
    i = 1;
    while J > Tolerance && i <= maxits
        Function = (1 / (M^2)) * ((f1 * (1 + (f2 * (M^2))))^(f3)) - (AR^2);
        dFunction = (f3 / (M^2)) * ((f1 * (1 + (f2 * (M^2))))^(f3 - 1)) * 2 * f2 * f1 * M + (-2 / (M^3)) * ((f1 * (1 + f2 * (M^2)))^(f3));
        y = M - (Function / dFunction);
        M = y;
        Function = (1 / (M^2)) * ((f1 * (1 + (f2 * (M^2))))^(f3)) - (AR^2);
        J = abs(Function);
        i = i + 1;
    end
end