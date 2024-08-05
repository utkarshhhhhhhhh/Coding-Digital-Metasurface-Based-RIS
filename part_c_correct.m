clc;
clear all;
close all;

%% Parameters

% Unit cell 
u_cell = 1.8e-3; 

% Reflectarray dimension
Array_dim = 54e-3 - u_cell; 

% Operating frequency (Hz)
frequency = 78.5e9;  

% Speed of light (m/s)
c = 3e8; 

% Periodicity of the reflectarray 
periodicity = 1.8e-3; 

% Feed point position 
z_feed = 80e-3;          

%% Derived variables

% Wavelength (m)
lambda = c / frequency;

% Wavenumber 
k0 = 2 * pi / lambda; 

% Target elevation angles (radians)
theta_targets = deg2rad([30, -30]);

% Target azimuth angle (radians)
phi_target = deg2rad(0);   

xf = -Array_dim/2 : u_cell : Array_dim/2;  % X positions across aperture (m)
yf = -Array_dim/2 : u_cell : Array_dim/2;  % Y positions across aperture (m)
[Xi, Yi] = meshgrid(xf, yf);  % Create grid of Xi and Yi positions

%% Initialize far-field pattern
theta = -pi/2:0.01:pi/2;  
E_far_field = zeros(size(theta));  

for theta_target = theta_targets
    %% Calculate the phase distribution across the array

    % Distance from feed to each point on the aperture (m)
    Ri = sqrt(Xi.^2 + Yi.^2 + z_feed^2);  
 
    % Spatial phase delay (radians)
    phi_spd = 0;  

    % Calculate the progressive phase contribution
    phi_pp = -k0 .* (Xi .* sin(theta_target) .* cos(phi_target) + Yi .* sin(theta_target) .* sin(phi_target));

    % Total phase distribution
    phi = phi_spd - phi_pp;  % Total phase (radians)

    % Convert phase to degrees for plotting
    phi_deg = rad2deg(phi);
    phi_deg = wrapTo360(phi_deg);

    %% Discretize phase to 1-bit (0° or 180°)

    phi_1bit = zeros(size(phi_deg));  % Initialize phi_1bit with zeros (same size as phi_deg)

    phi_1bit(phi_deg >= 90 & phi_deg < 270) = 180;  % Set 180° where phi_deg is in [90°, 270°)
    phi_1bit(phi_deg < 90 | phi_deg >= 270) = 0;    % Set 0° where phi_deg is not in [90°, 270°)

    %% Define parameters for far field calculation
    M = length(xf);  % Number of elements in the X dimension
    N = length(yf);  % Number of elements in the Y dimension

    % Convert phase to complex phasor for each element
    Gamma = 1 .* exp(1i * deg2rad(phi_1bit));
    
    % Initialize the resulting field E as a matrix of zeros
    E = zeros(M, N);

    % Calculate the field E
    for m = 1:M
        for n = 1:N
            E(m, n) = Gamma(m, n) * exp(-1j * k0 * ((m-1) * u_cell * sin(theta_target) * cos(phi_target) + (n-1) * u_cell * sin(theta_target) * sin(phi_target)));
        end
    end

    % Calculate and accumulate the far-field pattern
    for t = 1:length(theta)
        E_far_field(t) = E_far_field(t) + sum(sum(E .* exp(-1j * k0 * (Xi * sin(theta(t)) * cos(phi_target) + Yi * sin(theta(t)) * sin(phi_target)))));
    end
end

% Calculate the modulus of E(m,n)
disp(E);

% Normalize the far-field pattern
E_far_field_norm = abs(E_far_field) / max(abs(E_far_field));
E_far_field_dB = 20 * log10(E_far_field_norm);

%% Plot the normalized far-field pattern
theta_deg = rad2deg(theta);  
figure;
plot(theta_deg, E_far_field_dB, 'LineWidth', 2);
xlabel('\theta (degrees)');
ylabel('Normalized Pattern (dB)');
title('Normalized Far-Field Pattern for \theta_{target} = ±30 degrees');
grid on;
xlim([-90 90]);  
ylim([-30 0]);  
hold on;
