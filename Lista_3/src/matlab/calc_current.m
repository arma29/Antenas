function [current]=calc_current(lambda, length2,wire_radius,zprime,feed_voltage)
%function [current]=calc_current(length,radius,zprime,feed_voltage)
% lambda = wavelength
% length = length (in terms of wavelength) of the dipole
% radius = radius (in terms of wavelength) of the dipole
% zprime = zprime=linspace(-ant_length/2,ant_length/2,current_steps);
% feed_voltage = peak voltage at antenna feed.
%
% Solve for I in Z*I=V
% Where I is the unknown current distribution
% Z is a square matrix filled using Pocklingtons Integral Equation
% V is the user defined voltage (Becuase the antenna is center fed,
% only the middle entry has a value all other entries are zero.)

mu0 = 4*pi*10^-7; % permeability (in newtons/amp^2)
epsilon0 = 8.854187817*10^-12; % permittivity (in Farads/meter )
eta = sqrt(mu0/epsilon0); % intrinsic impedance (in ohms)
beta = (2*pi)/lambda; % wave number (in radians/meter)
current_steps = length(zprime); % of points for which to find a current
%
%1 2 3 4 5 6 7 8 9 A B C Point
% *-----*-----*-----*-----*-----*=====*-----*-----*-----*-----*-----*
%1 2 3 4 5 6 7 8 9 A B Segment
% ^
% gap
% The points of zprime are at the ends of each segment.
% zprime=linspace(-length2/2,length2/2,current_steps);
% The current is calculated at the center of each segment.
% This creates an index mismatch.
ie_steps = 199;
dz = zprime(2) - zprime(1);
ddz = dz/ie_steps;
% Caclulate a new vector with points at the center of each segment
%zcenter = zprime+.5*dz;
%zcenter = zcenter(1:length(zcenter)-1); %chop off the last element
%zdist = zcenter-zcenter(1); % distance from the end segment to every other segment
%rr = sqrt(wire_radius^2 + zdist.^2);

ZM = length2/2-.5*dz;
jj=1:ie_steps;
z = -.5*dz + ddz*jj;
ZN = 1:current_steps;
for ii=1:current_steps-1
    % crt = 0;
    ZN(ii) = length2/2 -(ii-.5)*dz;
    r = sqrt(wire_radius^2 + (ZN(ii)-ZM+z).^2); %distance from segment ii to skin of jj
    crt = pocklington(r,wire_radius, lambda)*ddz;
    p(ii) = sum(crt);
end
%
% Mark This code works without loops, but introduces a numerical "wobble" that throws-off
% the final directivity calculation by .02. I am sticking with the loop.--Rob
%keyboard
%P=p;
%ii=1:(current_steps-1);
%ZM = ZM+z;
%ZN=length2/2-(ii-.5)*dz;
%ZZ = ones(size(ii));
%ZZN = ZZ'*ZN;
%ZZM = ZM'*ZZ;
%r = sqrt(wire_radius^2 + (ZZN-ZZM).^2);
%crt = pocklington(r,wire_radius, lambda)*ddz;
%p = sum(crt);
%figure; plot(p-P); title('wobble P-p');

% Because Z is a toeplitz matrix we can fill in the entire matrix now that we have the
% first row. The general pattern is:
%
% |A B C D E F G . . . Z|
% |B A B C D E F . . . Y|
% |c B A B C D E . . . X|
% |. . .|
% |. . .|
% |. . .|
% |Z Y X W V U T . . . A|
znm=toeplitz(p);
middle = floor(current_steps/2);
G = zeros(size(p));
G(middle) = feed_voltage;
current=G/znm;

L = length(current)+1;
current(L)=1; % add the last point back to the array
temp = current(1:middle);
current(middle+1:L)=flipdim(temp,2);
