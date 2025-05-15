clear all
close all
clc

%% TENSIONED CABLE RESPONSE TO ASSIGNED INITIAL POSITION WITH ZERO VELOCITY
%% Definition of the mechanical properties of the system

load('init_pos.mat');
T=50000;                % tension [N]
m=14;                   % mass per unit length [kg/m]
c=sqrt(T/m);            % propagation velocity [m/s] 
L = dx*(length(x)-1);   % cable length [m]  

%% Setting the frequency range
fmax=5;                 %[Hz]
f=linspace(0,fmax,10^6);
omega=2*pi*f;           %[rad/s]

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

% Imposing that the determinant is null
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
    Ci_hat=[1; -Ei_hat\Hi_hat];
    
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
axis tight 

%% Visualization of the animation of one mode

figure(30), hold on, grid on, box on
title(['Mode ',num2str(mode)])
ylim([-5 5])
plot(x,phi(mode,:),':k','LineWidth',2)
h1=plot(x,zeros(size(x)),'LineWidth',2);
xlabel('Cable length [m]')
ylabel('Mode shape []')
axis tight

for t=linspace(0,2/f(i_nat(mode)),200)
    if ishandle(h1)
        w1=real(phi(mode,:)*exp(1i*omega(i_nat(mode))*t));
        h1.YData=w1;
        pause(.03)
    else
        return
    end
end

%% Initial position

% Plot initial position
figure;
hold on; box on; grid on;
ylim([-3 3])
h1 = plot(x,init_pos,'DisplayName','Initial Position')
legend

%% Response computation
% to be continued...
q = linspace(0,L,1001);
figure(100), hold on, grid on, box on
xlabel('Cable length [m]')
ylabel('Mode shape []')
ylim([-0.3,0.3])
h2 = zeros(16,1001);

for n = 1:16
   shape(n,:) = sin(n*pi/L*x).*init_pos;
   C2(n,1) =(2/L)*trapz(q,shape(n,:)); 
   h1(n,1)=plot(x,sin(n*pi/L*x)*C2(n,1),'LineWidth',2);
   h2(n,:) = sin(n*pi/L*x)*C2(n,1);
end

starting_point = sum(h2,1);

% questa sezione commentata contiene l'animazione delle standing waves che compongono la risposta del cavo
% for t=linspace(0,1,200) 
%     for m = 1:16
%       if ishandle(h1)
%         w1 = sin(m*pi/L*x)*C2(m,1)*cos(i_nat(1,m)*t);
%         h1(m,1).YData=w1;
%         pause(.03)
%       else
%         return
%       end
%     end
% end

figure(200), grid on, box on, hold on
xlabel('Cable length [m]')
ylabel('Time Response')
h3 = plot(x,starting_point,"LineWidth", 2,"Color", 'r');
plot(x, init_pos, 'LineStyle','--','Color', 'b')
ylim([-1,1])

for t=linspace(0,10,200)
    for k = 1:16
      if ishandle(h3)
        w(k,:) = sin(k*pi/L*x)*C2(k,1)*cos(omega(i_nat(1,k))*t);
      else
        return
      end
    end
    h3.YData = sum(w,1);
    pause(.10)

end

