%2D Array for calculating phase
clc;
clear all;
close all;

%% Parameters
 
% unit cell 
u_cell = 10e-3; 
% reflectarray dimension
Array_dim = 100e-3-u_cell;
% Operating frequency (Hz)
frequency = 1e9;  
% Speed of light (m/s)
c = 3e8; 
% Periodicity of the reflectarray 
periodicity = 10e-3; 
% Feed point position 
z_feed = 80e-3;          

%% Derived variables

% Wavelength (m)
lambda = c / frequency;
% Wavenumber 
k = 2 * pi / lambda; 
% Target elevation angle (radians)
theta_target = deg2rad(30); 
% Target azimuth angle (radians)
phi_target = deg2rad(0);     
xf =-Array_dim/2 : u_cell : Array_dim/2;  % X positions across aperture (m)
yf =-Array_dim/2 : u_cell : Array_dim/2;  % Y positions across aperture (m)
[Xi, Yi] = meshgrid(xf, yf);  % Create grid of Xi and Yi positions

%% Calculate the phase distribution across the 

% Distance from feed to each point on the aperture (m)
Ri = sqrt(Xi.^2 + Yi.^2 + z_feed^2);  
% Spatial phase delay (radians)
phi_spd = 0;  
% Calculate the progressive phase contribution
phi_pp = -k * (Xi * sin(theta_target) * cos(phi_target) + Yi * sin(theta_target) * sin(phi_target));
% Total phase distribution
phi = phi_spd - phi_pp;  % Total phase (radians)
% Convert phase to degrees for plotting
phi_deg = rad2deg(phi);


%% Pcolor plot of the phase distribution

figure(1);
pcolor(xf*1e3, yf*1e3,phi_deg);  % Plot in millimeters
shading flat;
colorbar;
title('Phase Distribution (Pcolor Plot)');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Phase (degrees)');

%% surf plot of the phase distribution

figure(2);
surf(xf*1e3, yf*1e3, phi_deg);  % Plot in millimeters
colorbar;
title('Phase Distribution (Surf Plot)');
xlabel('X Position (mm)');
ylabel('Y Position (mm)');
zlabel('Phase (degrees)');