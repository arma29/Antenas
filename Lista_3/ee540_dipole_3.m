function ee540_dipole_3
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
lambda = input('Wave Length (default=500): ');
if isempty(lambda)
lambda = 500; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %default value for lambda = 500
                                                                                                            end
                                                                                                            multiplier = input('\n\nAntenna Length as a fraction of Wavelength\n (.5 for 1/2
wavelength)\n (.005 for an infintesimal dipole)\n (1.25 to see some side
lobes)\n\nEnter Length Fraction (default=.5): ');
if isempty(multiplier)
multiplier = .5; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %default value for mutliplier = 0.5
                                                                                                           end
                                                                                                           ant_length = multiplier*lambda; %calculate the length of the antenna
multiplier2 = input('\n\nAntenna Radius as a fraction of Wavelength\n (.005 for thin
wire)\n 1/100 of dipole lenght is recommended \n\nEnter Radius Fraction (default=.005):
');
if isempty(multiplier2)
multiplier2 = .005; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %default value for mutliplier2 = 0.005
                                                                                                               end
                                                                                                               wire_radius = multiplier2*lambda; %calculate the diameter of the antenna
current_steps = input('\nNumber of constant current segments to use in approximating\n
the current (must be an odd number: default=199): ');
if isempty(current_steps)
current_steps=199; % % % % % % % % % % % % % % % % % % % % % % % % % % % % %default value for current_steps = 199
                                                                                                              end
                                                                                                              if (mod(current_steps,2)==1)
        current_steps = current_steps + 1; % ensure an odd number of segments & even number of
points
end
feed_voltage = input('\nPeak Voltage at Feed (default=.01*j): ');
if isempty(feed_voltage)
feed_voltage = .01*j; % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %default value for feed_voltage
end





%k is the wave number (called beta in some contexts).
k=2*pi/lambda;

%zprime holds the z-locations for the end of each current segment.
zprime=linspace(-ant_length/2,ant_length/2,current_steps);

%current holds an approximation of the current along each segment of the antenna.
%this code uses equation
current=calc_current3(lambda, ant_length, wire_radius, zprime, feed_voltage);
%plot the current distribution along the antenna
figure; plot(abs(current),zprime); title('Current Distribution (Magnitude)');
xlabel('Current (amps)'); ylabel('Z Location (m)');
%plot the current distribution along the antenna
figure; plot(angle(current),zprime); title('Current Distribution (Phase)');
xlabel('Current (amps)'); ylabel('Z Location (m)');

prad = 0;
T=[0:2*pi/360:2*pi];
P=[0:2*pi/360:2*pi];
F = dipoleint2(T, P, lambda, current, zprime);
%return
        %Calculate the total power radiated by "integrating" all values of the
        % radiation intensity times sin(theta)*dTheta*dphi
        % from theta = 0 to theta = pi;
% from phi=0 to 2*pi;
semicircle = sin(T(1 : 181)); %semicircle is a vector of sin(theta) from 0 to pi
prad1 = F(181,1 : 181); %at this point, prad1 sweeps through the ? ? Field from theta = 0 : pi
                                                                                        prad1 = prad1.*(pi/180).*semicircle; %now prad is F*sin(theta)*dtheta
prad = sum(prad1)*2*pi; %prad is the "area" under prad1 rotated 2*pi around phi.
%same prad calculation using loops for the double integration.
                             %prad=0;
%for (i=1 : 361)
        % for (j=1 : 181)
                % prad = prad + F(i,j)*sin(T(j))*(pi/180)^2;
% end
%end
M=4*pi*F./prad;
M1=M (91,:); %pluck out the graph that sweeps theta while phi=90
                                                               M2=M (:,91); %pluck out the graph that sweeps phi while theta=90
                                                                                                                              figure; polar (T,M1); title ('Elevation Plane');
figure; polar (P,M2); title ('Azmuthal Plane');
D0 = max (max (M));
Rrad = max (current);
Rrad = abs (2*prad/Rrad^2);
D0
prad
Rrad
        figure; mesh (M);
xlabel ('Elevation');
ylabel ('Azimuth');
zlabel ('Intensity (U)');
title (ant_length/lambda);
