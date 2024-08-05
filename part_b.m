clc;
clear;
close all;

%% Parameters
u_cell = 1.8e-3;      % Unit cell size (m)
Array_dim = 54e-3-u_cell;  % Reflectarray dimension (m)
frequency = 78.5e9;     % Operating frequency (Hz)
c = 3e8;             % Speed of light (m/s)
% periodicity = 1e-3; % Periodicity of the reflectarray (m)
z_feed = 80e-3;      % Feed point position (m)


%% Derived variables

lambda = c / frequency;  % Wavelength (m)
k = 2 * pi / lambda;     % Wavenumber

theta_target = deg2rad(30);  % Target elevation angle (radians)
phi_target = deg2rad(0);     % Target azimuth angle (radians)

xf = -Array_dim/2 : u_cell : Array_dim/2;  % X positions across aperture (m)
yf = -Array_dim/2 : u_cell : Array_dim/2;  % Y positions across aperture (m)
[Xi, Yi] = meshgrid(xf, yf);  % Create grid of Xi and Yi positions

%% Calculate the phase distribution across the array

Ri = sqrt(Xi.^2 + Yi.^2 + z_feed^2);  % Distance from feed to each point on the aperture (m)

phi_spd = 0;  % Spatial phase delay (radians)

phi_pp = -k * (Xi * sin(theta_target) * cos(phi_target) + Yi * sin(theta_target) * sin(phi_target));  % Progressive phase contribution

phi = phi_spd - phi_pp;  % Total phase (radians)
phi_deg = rad2deg(phi);  % Convert phase to degrees for plotting
phi_deg = wrapTo360(phi_deg);

%% Discretize phase to 1-bit (0° or 180°)

phi_1bit = zeros(size(phi_deg));  % Initialize phi_1bit with zeros (same size as phi_deg)

phi_1bit(phi_deg >= 90 & phi_deg < 270) = 180;  % Set 180° where phi_deg is in [90°, 270°)
phi_1bit(phi_deg < 90 | phi_deg >= 270) = 0;    % Set 0°where phi_deg is not in [90°, 270°)

%% Pcolor plot of the phase distribution

figure(1);
pcolor(xf*1e3, yf*1e3, phi_1bit);  % Plot in millimeters
 
shading flat;
colorbar;

clim([0 180]);  % Set color axis limits
title('1-bit Phase Distribution (Pcolor Plot)');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Phase (degrees)');

%% Surf plot of the phase distribution

figure(2);
surf(xf*1e3, yf*1e3, phi_1bit);  % Plot in millimeters

colorbar;

clim([0 180]);  % Set color axis limits
title('1-bit Phase Distribution (Surf Plot)');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Phase (degrees)');
