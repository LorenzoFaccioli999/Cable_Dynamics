clear all
close all
clc

%% TENSIONED CABLE 
%% Definition of the mechanical properties of the system

T=500000;               % tension [N]
m=14;                   % mass per unit length [kg/m]
c=sqrt(T/m);            % propagation velocity [m/s] 
L = 100;                % cable length [m]  

%% Setting the frequency range
fmax=5;                 %[Hz]
f=linspace(0,fmax,10^6);
omega=2*pi*f;           %[rad/s]

%% Setting the space domain
x=linspace(0,L,10^3);

%% Building the matrix of the coefficients from the BCs

H=@(omega) [  0             1    ;
    sin(omega/c*L)            cos(omega/c*L)];

%% Plotting the determinat as a function of omega and natural frequencies

for i=1:length(omega)
    dets(i)=det(H(omega(i)));
end

figure(10), box on
semilogy(f,abs(dets))
hold on, grid on, xlabel('f [Hz]')
title(['Natural frequencies in the range 0 - ', num2str(fmax), ' Hz'])

% Enforcing that the determinant is null
i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end

plot(f(i_nat),abs(dets(i_nat)),'or','Linewidth',1)

%% Solving the reduced system

for i_mode=1:length(i_nat)
    omega_i=omega(i_nat(i_mode));
    Hi=H(omega_i);
    Hi_hat=Hi(2,2);
    Ei_hat=Hi(2,1);
    Ci_hat=[1; -Hi_hat\Ei_hat];
    
    C_hat(:,i_mode)=Ci_hat;
end

%% Mode shapes computation

for i_mode=1:length(i_nat)
    omega_i=omega(i_nat(i_mode));
    phi(i_mode,:)= C_hat(1,i_mode)*sin(omega_i/c*x) + C_hat(2,i_mode)*cos(omega_i/c*x);
end

prompt2={'Enter the mode shape to display:'};
answer2=inputdlg(prompt2);
mode=str2double(answer2);

figure(20)
plot(x,phi(mode,:),'-k','Linewidth',2)
ylim([-5 5])
grid on
xlabel('Cable length [m]')
ylabel('Mode shape []')
title(['Mode ',num2str(mode)])

%% Visualization of the animation of one mode

figure(30), hold on, grid on, box on
tit = sprintf('Mode %d - f%d = %.2f Hz', mode, mode, f(i_nat(mode)));
title(tit)
ylim([-5 5])
plot(x,phi(mode,:),':k','LineWidth',2)
h1=plot(x,zeros(size(x)),'LineWidth',2);
xlabel('Cable length [m]')
ylabel('Mode shape []')

for t=linspace(0,2/f(i_nat(mode)),200)
    if ishandle(h1)
        w1=real(phi(mode,:)*exp(1i*omega(i_nat(mode))*t));
        h1.YData=w1;
        pause(.03)
    else
        return
    end
end

%% Assigned position

wo = sin(2*pi/L*x);

% Plot initial position
figure;
hold on; box on; grid on;
ylim([-3 3])
h2 = plot(x,wo,'DisplayName','Progresive wave', 'LineWidth', 2)
h3 = plot(x,wo,'DisplayName','Regressive wave', 'LineWidth', 2)
h4 = plot(x,wo,'DisplayName','Second mode shape', 'LineWidth', 6)
legend

%% Response computation with travelling wave solution
% to be continued...
for t=linspace(0,2/f(i_nat(mode)),200)
    if ishandle(h2)
        w2=0.5*sin(2*pi/L*(x-c*t));
        w3= 0.5*sin(2*pi/L*(x+c*t));
        w4 = w2+w3;
        h2.YData= w2;
        h3.YData = w3;
        h4.YData= w4;
        pause(.03)
    else
        return
    end
end



