%Aaron Rosen - Fields & Waves II - Project 4 - FD1D_13.c to MATLAB

%Modify your program fd1d_1.3.c to simulate the sinusoidal source (see
%fd1d_1.4.c)

%Keep increasing your incident frequency from 700MHz upward at intervals of
%300MHz. What happens?

%A type of propagating wave function that is of great interest in areas
%such as optics is the "wave packet," which is a sinusoidal function in a
%Gaussian envelope. Modify your program to simulate a wave packet.

KE = 200;
Ex = zeros(KE, 1);
Hy = zeros (KE, 1);
CB = zeros(KE, 1);

ex_low_m1 = 0;
ex_low_m2 = 0;
ex_high_m1 = 0;
ex_high_m2 = 0;

spread = 121
t0 = 40.0;
T = 0;
NSTEPS = 1;
kc = KE/2;
ddx =0.01;
dt = ddx/(2*3e8);


for k = 2:KE
    CB(k) = 0.5;
end

dielectric = "Dielectric starts @ ---> ";
kstart = input(dielectric);
disp(kstart);
eps = "Input Epsilon ---> ";
epsilon = input(eps);
disp(epsilon);
frequency = "Input frequency (MHz) ---> ";
freq = input(frequency)*1e6;

for k =kstart:KE
    CB(k) = 0.5/epsilon;
end


while (NSTEPS > 0)
    n = 0;
    prompt = "NSTEPS ---> ";
    NSTEPS = input(prompt);
    disp(NSTEPS);
    
    for n=1:NSTEPS
        T = T+1;

        %Calculuate E-field
        for k =2:KE
            Ex(k) = Ex(k) + CB(k)*(Hy(k-1) - Hy(k));
        end

               %Define Boundary Conditions
        Ex(1) = ex_low_m2;
        ex_low_m2 = ex_low_m1;
        ex_low_m1 = Ex(2);
        Ex(KE-1) = ex_high_m2;
        ex_high_m2 = ex_high_m1;
        ex_high_m1 = Ex(KE-2);


        %Add Sinusoidal Gaussian Pulse @ low end
        pulse = sin(2*pi*freq*dt*T);
        Ex(6) = Ex(6) + pulse;
        disp(t0-T);
        disp(Ex(6));

    

    for k = 1:KE-1
        Hy(k) = Hy(k) + 0.5* (Ex(k) - Ex(k+1));
    end

    figure(1);
    subplot(2,1,1)
    plot(Ex, 'LineWidth', 2)
    ylabel('Ex')
    ylim([-1 1])

    subplot(2,1,2)
    plot(Hy, LineWidth=2);
    ylabel('Hy')
    ylim([-5 5])
    for k = 1:KE
        disp(Ex(k));
    end

    end
    cont = "Would you like to continue?";
    if (input(cont) == 0)
        clear all; close all; clc;
    else if (input(cont) == 1)
            continue;
    end
    end
    

end

