%% CABLE WITH DIFFERENT MASS PER UNIT LENGTH ON THE TWO HALVES - Free and forced motion in freq. range 0-4.5Hz

%%
clear all
close all
clc

%% Input geometrical data
% We consider a 100 metres long steel cable, 
% with different mass per unit length on the two halves.
% 
% We define the physical properties of the beam:

L = 50;               % m       length of each subspan
m2 = 14;              % Kg/m    mass per unit length
m1 = 3*m2;            % Kg/m    mass per unit length
T = 100e3*9.81;       % N       tension force

c1 = sqrt(T/m1);
c2 = sqrt(T/m2);
F0 = 10000;

%% Setting the frequency range
% We investigate the frequency range from 0 to 4.5 Hz
fmax = 4.5;
f=linspace(0,4.5,2000); %ho ristretto la lunghezza del vettore perchè 
                        %impiegava troppo tempo altrimenti
omega=2*pi*f;
x = linspace(0,L,length(f));
x_tot = linspace(0,2*L,length(f)*2);

%% Setting the boundary conditions
% Suggestion: consider the reference systems as the one suggested in the
% lectures... (1st: from left to right, 2nd: for right to left)

%% Building the H matrix
H=@(omega) [  sin((omega/c1)*L)             -sin((omega/c2)*L)    ;
           (T*omega/c1)*cos((omega/c1)*L)  (T*omega/c2)*cos((omega/c2)*L)];

%% Plotting the determinat as a function of omega
for i=1:length(omega)
    dets(i) = det(H(omega(i)));
end

figure(10), box on
semilogy(f,abs(dets))
hold on, grid on, xlabel('f [Hz]')
title(['Natural frequencies in the range 0 - ', num2str(fmax), ' Hz'])


%% Imposing that the determinant is null
i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end

plot(f(i_nat),abs(dets(i_nat)),'or','Linewidth',1)

%% Solving the reduced system
% Now we know the value for which the system is singular (i.e. admits 
% non-trivial solutions), we can find the modal shapes solving the reduced system
F_vect = [0 ; F0];

for jj=2:length(omega)
     H_j = H(omega(jj));
     A1 = F0/(H_j(2,1)-H_j(1,1)*H_j(2,2)/H_j(1,2));
     A2 = -H_j(1,1)*A1/H_j(1,2);
%Cramer, meno efficiente in questo caso
%     A1_matrix = [F_vect, [H_j(1,2); H_j(2,2)]];
%     A2_matrix = [[H_j(1,1); H_j(2,1)], F_vect];
%     A1 = det(A1_matrix)/det(H_j);
%     A2 = det(A2_matrix)/det(H_j);

  
    C(:,jj)=[A1,A2]; %vettore dei coefficienti
end

%% Mode shapes computation
for mode=2:length(omega)
    omega_i=omega(mode);
    phi1(mode,:)= C(1,mode)*sin(omega_i/c1*x);
    phi2(mode,:) =  C(2,mode)*sin(omega_i/c2*x);
    phi_12(mode,:) = [phi1(mode,:) , fliplr(phi2(mode,:))];
end

prompt2={'Enter the mode shape to display:'};
answer2=inputdlg(prompt2);
mode=str2double(answer2);

figure(20)
plot(x_tot,phi_12(i_nat(mode),:),'-k','Linewidth',2)
% ylim([-100 100])
grid on
xlabel('Cable length [m]')
ylabel('Mode shape []')
title(['Mode ',num2str(mode)])


%% Visualization of the animation
figure(30), hold on, grid on, box on
title(['Mode ',num2str(mode)])
% ylim([-5 5])
plot(x_tot,phi_12(i_nat(mode),:),':k','LineWidth',2)
h1=plot(x_tot,zeros(size(x_tot)),'LineWidth',2);
xlabel('Cable length [m]')
ylabel('Mode shape []')

for t=linspace(0,2/f(i_nat(mode)),200)
    if ishandle(h1)
        w1=real(phi_12(i_nat(mode),:)*exp(1i*omega(i_nat(mode))*t));
        h1.YData=w1;
        pause(.03)
    else
        return
    end
end

%% FORCED MOTION: Solving the non-homogeneous system
% Plot of the values of the coefficients A1 and A2 as a function of Omega
figure(40), hold on, grid on, box on, axis tight
title('A1 A2')
xlabel('Omega')
plot(omega, C(1,:),'LineWidth',1,'Color','r', 'DisplayName', 'A1')
plot(omega, C(2,:),"LineWidth",1,"Color",'b', 'DisplayName', 'A2')
legend
ylim([-3 3])
%% Response computation for a selected value of Omega
prompt3={'Enter the desired value for omega:'};
answer3=inputdlg(prompt3);
omega_selected=str2double(answer3);
f_selected = omega_selected/(2*pi);
H_s = H(omega_selected);
A1_matrix_s = [F_vect, [H_s(1,2); H_s(2,2)]];
A2_matrix_s = [[H_s(1,1); H_s(2,1)], F_vect];
A1_s = det(A1_matrix_s)/det(H_s);
A2_s = det(A2_matrix_s)/det(H_s);
  
C_s=[A1_s;A2_s];

phi1_s= C_s(1,1)*sin(omega_selected/c1*x);
phi2_s =  C_s(2,1)*sin(omega_selected/c2*x);
phi_12_s = [phi1_s , fliplr(phi2_s)];

figure(50), hold on, grid on, box on, axis tight
title(['Response of the cable'])
% ylim([-0.3 0.3])
plot(x_tot,phi_12_s,':k','LineWidth',2)
h2=plot(x_tot,zeros(size(x_tot)),'LineWidth',2);
xlabel('Cable length [m]')
ylabel('Mode shape []')



%% Visualization of the animation
for t=linspace(0,6/f_selected,200)
    if ishandle(h2)
        w2=real(phi_12_s*exp(1i*omega_selected*t));
        h2.YData=w2;
        pause(.03)
    else
        return
    end
end


