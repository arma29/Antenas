
% author: UJWOL PALANCHOKE
%         JACOBS UNERSITY, BREMEN, GERMANY
% all righ reserved to author

% Impedance of coupled dipole antenna with respect to wavelength.The dipole antenna seperated by distance
% delta. The method used is Hallen's Integral. Here the second dipole is not fed by any
% voltage source


clear all;
lamda=1;        %arbitary wavelength (in meter): this can be varied accourding to the problem

l=0.5*lamda;    % length of dipole: this can be varied accourding to the problem
a=0.001*lamda;  % radius of dipole: this can be varied accourding to the problem
N=50;           % no.of segments i.e. the sampling points in antenna surface/for convergence: this can be varied accourding to the problem
dz=l/(2*(N-1/2));% incremental length  of l/2 upper part of antenna
z=dz/2:dz:dz*(N); % stores the invremental length in upper part of antenna
k=2*pi/lamda;       % wave number
etta=377;           % free space impedance
V_impressed=-1;      % impressed voltage on D1: this can be varied accourding to the problem

delta=0.1*lamda:0.1:5*lamda;  % distance between the dipoles

for d=1:length(delta)
%%%%%% calculating the Green's Function  on first dipole (with source)%%%%%

%%%%%%%%%%%%%%%%calculating Function's due to current in first and second dipole%%%%%%%%

for m=1:N %m is the obsevation point
    for n=1:N-1% n is the source points
        R1=sqrt(a^2+(z(m)-z(n))^2);     % R+ in lecture notes for G11
        R2=sqrt(a^2+(z(m)+z(n))^2);     % R- in lecture notes for G11
        R3=sqrt(delta(d)^2+(z(m)-z(n))^2); % R+ for G12
        R4=sqrt(delta(d)^2+(z(m)+z(n))^2); % R- for G12
        
%%%%%%%%%%%%% evaluating greens function G11 %%%%%%%%%%%%%%%%
        G1=exp(-j*k*R1)/(4*pi*R1);      % G11(R+) in lecture notes
        G2=exp(-j*k*R2)/(4*pi*R2);      % G11(R-) in lecture notes
        G_11(m,n)=(G1+G2)*dz;           % G11(R) observation in D1 source:1
%%%%%%%%%%%%% evaluating greens function G12 %%%%%%%%%%%%%%
         G3=exp(-j*k*R3)/(4*pi*R3);      % G12(R+) in lecture notes
         G4=exp(-j*k*R4)/(4*pi*R4);      % G12 (R-)in lecture notes
         G_12(m,n)=(G3+G4)*dz;           % G12(R) observation in D1 source:2
    end
    v_1(m)=j/(2*etta)*(sin(k*z(m)));      % Sine term due to impressed voltage in D1
end

%%%%%% calculating the Green's Function  on second  dipole (without source)%%%%%

%%%%%%%%%%%%%%%%calculating Function's due to current in first and second dipole%%%%%%%%

for m=1:N %m is the obsevation point
    for n=1:N-1% n is the source points
        R5=sqrt(a^2+(z(m)-z(n))^2);     % R+ in lecture notes for G22
        R6=sqrt(a^2+(z(m)+z(n))^2);     % R- in lecture notes for G22
        R7=sqrt(delta(d)^2+(z(m)-z(n))^2); % R+ for G21
        R8=sqrt(delta(d)^2+(z(m)+z(n))^2); % R- for G21
        
%%%%%%%%%%%%% evaluating greens function G22 %%%%%%%%%%%%%%
        G5=exp(-j*k*R5)/(4*pi*R5);      % G22(R+) in lecture notes
        G6=exp(-j*k*R6)/(4*pi*R6);      % G22(R-) in lecture notes
        G_22(m,n)=(G5+G6)*dz;           % G22(R) observation in D2 source:2
%%%%%%%%%%%%% evaluating greens function G12 %%%%%%%%%%%%%%
         G7=exp(-j*k*R7)/(4*pi*R7);      % G21(R+) in lecture notes
         G8=exp(-j*k*R8)/(4*pi*R8);      % G21 (R-)in lecture notes
         G_21(m,n)=(G7+G8)*dz;           % G21(R) observation in D2 source:1
    end
    v_2(m)=0;      % Sine term due to impressed voltage in D2 (no voltage applied in D2)
end

%%%%%%%%%%%% concatenating G11 and G12%%%%%%%%%%%
G_D1=horzcat(G_11,G_12);
s=size(G_D1);
N_t1=s(2);
G_D1(:,N_t1+1)=[-cos(k*z)];        % N element holds the cosine terms for D1
G_D1(:,N_t1+2)=[0];              % N+1 element holds the cosine terms for D2 which is zero for D1


%%%%%%%%%%%% concatenating G21 and G22%%%%%%%%%%%
G_D2=horzcat(G_21,G_22);
s_1=size(G_D2);
N_t2=s_1(2);
G_D2(:,N_t2+1)=[0];           % N element holds the cosine terms from D1 which is zero for D1
G_D2(:,N_t2+2)=[-cos(k*z)];   % N+1 element holds the cosine terms from D2 

%%%%%%%%%%%%%concatenating Sine terms i.e. Vm%%%%%%%%%%%%%
Vm=horzcat(v_1,v_2);

%%% Defining the matrix required to evaluate the current%%%%%%%%
G_eval=vertcat(G_D1,G_D2);   
%%%% evaluate current%%%

I=inv(G_eval)*Vm.';     % this will give a column vector of current

%%%%%%%%%%% for question number 4 (LAB part)%%%%%%%%%%%%%%%%%%%%
s_2= size(I);
N_t3=s_2(1);               % this will give the number of elements in column
N1=(N_t3-2)/2;
I1=I(1);            % first element gives the current in the centre of D1
I2=I(N1+1);         % this gives the current in the centre of D2

Z_11(d)=V_impressed/I1;
Z_12(d)=V_impressed/I2;
end
figure(1)
plot(z,I);
xlabel('teste');
ylabel('teste2');
grid on

figure(2)

plot(delta,real(Z_11));
hold on
plot(delta,imag(Z_11),'r')
legend('real-Z11','Imaginary-Z11')
title('real and imaginary part of impedance')
xlabel('delta');
ylabel('impedance')
grid on

figure(3)

plot (delta, real(Z_12),'b:<')
hold on
plot (delta,imag(Z_12),'r-->')
legend('real-Z12','Imaginary-Z12')
title('real and imaginary part of impedance')
xlabel('delta');
ylabel('impedance')
grid on


    
        
        


