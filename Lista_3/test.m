clc;
clear all;
close all;

lambda = 1;
L = lambda/2;
a = lambda*10e-4;
N = 6;
V0 = 1;

k=2*pi/lambda;

Z=linspace(-L/2,L/2,N+1);
V = zeros(1,N)' %o ' da transpose, virando um vetor coluna
MtzZ = zeros([N])



for i=1:N-1
    V(i) = 1;
    if(i==1)
        V(i) = 0;
    end
end