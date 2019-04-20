function ee540_dipole_3
clc;
clear all;
close all;
% ee540_dipole3--
% Calculate the Directivity (D0), Radiation Resistance (Rrad), and
% Radiated Power (prad) for a given wavelength (lambda), dipole length,
% and dipole diameter.
%
% This calculation is for a dipole antenna of any length,
% centered at the origin and running along the z-axis.
%
% Prompt the User for lambda, antenna length (in terms of lambda),
% antenna radius (in terms of lambda), the number
% of current steps, and feed voltage.
lambda = 1;
ant_length = lambda/2;
wire_radius = lambda*10e-4;
current_steps = 60;
feed_voltage = 1;

%k is the wave number (called beta in some contexts).
k=2*pi/lambda;

%zprime holds the z-locations for the end of each current segment.
zprime=linspace(-ant_length/2,ant_length/2,current_steps);

%current holds an approximation of the current along each segment of the antenna.
%this code uses equation
current=calc_current(lambda, ant_length, wire_radius, zprime, feed_voltage);
%plot the current distribution along the antenna
figure; plot(abs(current),zprime); title('Current Distribution (Magnitude)');
xlabel('Current (amps)'); ylabel('Z Location (m)');
%plot the current distribution along the antenna
figure; plot(angle(current),zprime); title('Current Distribution (Phase)');
xlabel('Current (amps)'); ylabel('Z Location (m)');


