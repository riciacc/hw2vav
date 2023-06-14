clear all, close all, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          HOMEWORK 2                                               %
%                   Iaccarino Riccardo, Dalla Mora Enrico, Scarano Rocco                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data & Matrices    updated data
M1 = 5;
J1 = 2.5;
R1 = 1;

M2 = 1.25;
J2 = 0.16;
R2 = 0.5;

M3 = 10;

k1 = 1000;
c1 = 0.5;

k2 = 100;
c2 = 0.5;

k3 = 560;
c3 = 1;

k4 = 800;
c4 = 4;



%create the matrices

jac_M = [0  R2  1
         1  0  0
         0  0  1
         0  1  0
         0  0  1];

M = diag([M1 J1 M2 J2 M3]);
M_gen = jac_M'*M*jac_M;

jac_K = [0  -R2  -1
        -R1  2*R2  0
         R1  R2  0
         0  0  1];

K = diag([k1 k2 k3 k4]);
K_gen = jac_K'*K*jac_K;

jac_C = jac_K;
C = diag([c1 c2 c3 c4]);
C_gen = jac_C'*C*jac_C;

%init cond
x_3_0 = 0.1;
theta_1_0 = pi/12;
theta_2_0 = -pi/12;
x_dot_3_0 = 1;
theta_dot_1_0 = 0.5;
theta_dot_2_0 = 2;

%% UNDAMPED - Eigenfrequencies & eigenvectors

% eigenvalues eigenvectors problem
% D eigenvalues
% V eigenvectors
[V,D] = eig(inv(M_gen)*K_gen);
w_nat = sqrt(diag(D));
% normalization on the first independent variable
V_1 = [V]./[V(1,:)];


w_nat_ord= sort(w_nat);
w_nat_ord = [w_nat_ord(1); -w_nat_ord(1); w_nat_ord(2); -w_nat_ord(2); w_nat_ord(3); -w_nat_ord(3)];

V_ord = [V_1(:,3),  V_1(:,3), V_1(:,2), V_1(:,2), V_1(:,1), V_1(:,1)];

%% DAMPED - Eigenfrequencies & eigenvectors

% in [A] the submatrix [C] is full
A_damp = - [M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen]\[C_gen K_gen; -M_gen zeros(size(M_gen))];
[V_damp,D_damp] = eig(A_damp);
V_damp = V_damp(4:6,:); %consider only the displacements - get rid of the speed
lambda = diag(D_damp); %complex and conjugate
w_d_ord = [imag(lambda(1)); imag(lambda(3)); imag(lambda(5))];
[w_d_ord, index] = sort(w_d_ord);
w_d_ord = [w_d_ord(1); -w_d_ord(1); w_d_ord(2); -w_d_ord(2); w_d_ord(3); -w_d_ord(3)];
V_damp_ord = [V_damp(:,index(1)*2-1), V_damp(:,index(1)*2), V_damp(:,index(2)*2-1), V_damp(:,index(2)*2), V_damp(:,index(3)*2-1), V_damp(:,index(3)*2)];
V_damp_norm = [V_damp_ord]./[V_damp_ord(1,:)];

%% Rayleigh damping

% Using not ordered w_nat and w_d
a = real(lambda);
w_d = imag(lambda);
csi = abs(a./sqrt(w_d.^2+a.^2)); %equals |a/w_0|
h = [csi(1); csi(3); csi(5)];
A = [(1/2).*(1./w_nat), w_nat./2];
const_mod = A\h;
alpha = const_mod(1);
beta = const_mod(2);

% % Using ordered w_nat and w_d
% a = [real(lambda(index(1)*2-1)); real(lambda(index(1)*2)); real(lambda(index(2)*2-1)); real(lambda(index(2)*2)); real(lambda(index(3)*2-1)); real(lambda(index(3)*2))];
% csi = abs(a./sqrt(w_d_ord.^2+a.^2));
% h = [csi(1); csi(3); csi(5)];
% w_0= sort(w_nat);
% A = [(1/2).*(1./w_0), w_0./2];
% const_mod = A\h;
% alpha1 = const_mod(1);
% beta1 = const_mod(2);

% manually evaluated alpha and beta
% Sm =0;
% Skm =0;
% Scm =0;
% Sck =0;
% Sk = 0;
% for i = 1:9
%     Sm = Sm + M_gen(i)^2;
%     Skm = Skm + M_gen(i)*K_gen(i);
%     Scm = Scm + C_gen(i)*M_gen(i);
%     Sck = Sck + C_gen(i)*K_gen(i);
%     Sk = Sk + K_gen(i)^2;
% end

%const_mod = lsqr( [M_gen(:) , K_gen(:)] , C_gen(:)); %pseudo-inverse (least mean square error)


% beta_man = ((Sck/Sk)-((Scm*Skm)/(Sk*Sm)))/(1-((Skm^2)/(Sk*Sm)));
% 
% alpha_man = (Scm - beta_man*Skm)/Sm;

%C_gen_ray = alpha_man*M_gen + beta_man*K_gen;
C_gen_ray = alpha*M_gen + beta*K_gen;

%% Evaluation of damped frequencies considering the Rayleigh damping
A_damp_1 = - inv([M_gen zeros(size(M_gen)); zeros(size(M_gen)) M_gen])*[C_gen_ray K_gen; -M_gen zeros(size(M_gen))];
[V_damp_1,D_damp_1] = eig(A_damp_1);
V_damp_1 = V_damp_1(4:6,:); %consider only the displacements - get rid of the speed
lambda_1 = diag(D_damp_1); %complex and conjugate
w_d_ord_1 = [imag(lambda_1(1)); imag(lambda_1(3)); imag(lambda_1(5))];
[w_d_ord_1, index_1] = sort(w_d_ord_1);
lambda_1_ord = [lambda_1(index_1(1)*2-1); lambda_1(index_1(1)*2);
              lambda_1(index_1(2)*2-1); lambda_1(index_1(2)*2);
              lambda_1(index_1(3)*2-1); lambda_1(index_1(3)*2)];
w_d_ord_1 = [w_d_ord_1(1); -w_d_ord_1(1); w_d_ord_1(2); -w_d_ord_1(2); w_d_ord_1(3); -w_d_ord_1(3)];
V_damp_ord_1 = [V_damp_1(:,index_1(1)*2-1), V_damp_1(:,index_1(1)*2), V_damp_1(:,index_1(2)*2-1), V_damp_1(:,index_1(2)*2), V_damp_1(:,index_1(3)*2-1), V_damp_1(:,index_1(3)*2)];
V_damp_norm_1 = [V_damp_ord_1]./[V_damp_ord_1(1,:)];
V_damp_norm_1 =real(V_damp_norm_1);
%% 2.1 Free motion of the system

t = 0:0.01:20; %time axis

pos_t0 = [theta_1_0; theta_2_0; x_3_0]; %vector of the initial positions
vel_t0 = [theta_dot_1_0; theta_dot_2_0; x_dot_3_0]; %vector of the initial velocities

const = [V_damp_norm_1;lambda_1_ord.'.*V_damp_norm_1]\[pos_t0;vel_t0];

free_motion = zeros(3,length(t));

for ii = 1:length(t)
    free_motion(:,ii) = (const.'.*V_damp_norm_1)*exp(lambda_1_ord*t(ii));
end

%Plot of Free Damped Response

figure(1)
newcolors = [0.3 0.647 0.891
             0.95 0.625 0.298
             0.625 0.3 0.6];
        colororder(newcolors);
subplot(311)
plot(t,real(free_motion(1,:)), 'Color', newcolors(1,:), LineWidth=1.1);
xlabel('Time [s]');
ylabel({'Rotation \theta_1 [rad]'});
ylim([-0.3, 0.3]);
grid on
title('\theta_1')


subplot(312)
plot(t,real(free_motion(2,:)),'Color', newcolors(2,:), LineWidth=1.1);
xlabel('Time [s]');
ylabel({'Rotation \theta_2 [rad]'});
ylim([-0.3, 0.3]);
grid on
title('\theta_2')


subplot(313)
plot(t,real(free_motion(3,:)),'Color', newcolors(3,:), LineWidth=1.1);
xlabel('Time [s]');
ylabel({'Displacement x_3 [m]'});
grid on
title('x_3')


%% 2.2 Set initial condition for single mode contribution

pos_t0_1st_mode = [V_damp_norm_1(:,5)]; %vector of the initial positions
vel_t0_1st_mode = [0;0;0]; %vector of the initial velocities 
%vel_t0_1st_mode = [V_damp_norm_1(:,1)]; % with this set of velocities we
%obtain quite the same results, because theese number are consistent with
%the third mode of vibration
const_1st_mode = inv([V_damp_norm_1; lambda_1_ord.'.*V_damp_norm_1])* [pos_t0_1st_mode;vel_t0_1st_mode];

free_motion_1st_mode = zeros(3,length(t));

for ii = 1:length(t)
    free_motion_1st_mode(:,ii) = (const_1st_mode.'.*V_damp_norm_1)*exp(lambda_1_ord*t(ii));
end

%Plot of Free Damped Response

figure(2),
subplot 311
plot(t,real(free_motion_1st_mode(1,:)),'r', LineWidth=1.1);
xlabel('Time [s]');
ylabel({'displacement y_1','[m]'});
grid on
title('y_1')


subplot 312
plot(t,real(free_motion_1st_mode(2,:)),'b', LineWidth=1.1);
xlabel('Time [s]');
ylabel({'Rotation \theta_2','[rad]'});
grid on
title('\theta_2')


subplot 313
plot(t,real(free_motion_1st_mode(3,:)),'m', LineWidth=1.1);
xlabel('Time [s]');
ylabel({'Rotation \theta_3','[rad]'});
grid on
title('\theta_3')




%prove sulle initial conditions

% %% 2.2 Set initial condition for single mode contribution
% 
% pos_t0_1st_mode = [10;20;-30]; %vector of the initial positions
% vel_t0_1st_mode = [lambda_1_ord(5)*V_damp_norm_1(:,5)]; %vector of the initial velocities 
% %vel_t0_1st_mode = [V_damp_norm_1(:,1)]; % with this set of velocities we
% %obtain quite the same results, because theese number are consistent with
% %the third mode of vibration
% const_1st_mode = [V_damp_norm_1; lambda_1_ord.'.*V_damp_norm_1]\[pos_t0_1st_mode;vel_t0_1st_mode];
% 
% free_motion_1st_mode = zeros(3,length(t));
% 
% for ii = 1:length(t)
%     free_motion_1st_mode(:,ii) = (const_1st_mode.'.*V_damp_norm_1)*exp(lambda_1_ord*t(ii));
% end
% 
% %Plot of Free Damped Response
% 
% figure(40),
% subplot 311
% plot(t,real(free_motion_1st_mode(1,:)),'r', LineWidth=1.1);
% xlabel('Time [s]');
% ylabel({'displacement y_1','[m]'});
% grid on
% title('y_1')
% 
% 
% subplot 312
% plot(t,real(free_motion_1st_mode(2,:)),'b', LineWidth=1.1);
% xlabel('Time [s]');
% ylabel({'Rotation \theta_2','[rad]'});
% grid on
% title('\theta_2')
% 
% 
% subplot 313
% plot(t,real(free_motion_1st_mode(3,:)),'m', LineWidth=1.1);
% xlabel('Time [s]');
% ylabel({'Rotation \theta_3','[rad]'});
% grid on
% title('\theta_3')
% 
% sgtitle('Damped time responce with the only contribution of the third mode')

%% 3.a Forced Motion

w=0:0.01:30;

for ii=1:length(w)
    H = inv(-w(ii)^2*M_gen+1i*w(ii)*C_gen_ray+K_gen);
    H11(ii)=H(1,1);
    H12(ii)=H(1,2);
    H13(ii)=H(1,3);
    H21(ii)=H(2,1);
    H22(ii)=H(2,2);
    H23(ii)=H(2,3);
    H31(ii)=H(3,1);
    H32(ii)=H(3,2);
    H33(ii)=H(3,3);

end
y_ticks = [-180, -90, 0, 90, 180];
y_lim=[-5,5];

figure(31)
set(gcf,'position',[70, 70, 1400, 500])
subplot(2,3,1); plot(w,abs(H11),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_,_1'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,2); plot(w,abs(H12),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_,_2'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,3); plot(w,abs(H13),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_,_3'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,4); plot(w,angle(H11),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_1) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,5); plot(w,angle(H12),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_2) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,6); plot(w,angle(H13),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_3) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});

figure(32)
set(gcf,'position',[70, 70, 1400, 500])
subplot(2,3,1); plot(w,abs(H21),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [rad/N]'); title('H_2_,_1'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,2); plot(w,abs(H22),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [rad/N]'); title('H_2_,_2'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,3); plot(w,abs(H23),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_3| [rad/N]'); title('H_2_,_3'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,4); plot(w,angle(H21),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_1) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,5); plot(w,angle(H22),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_2) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,6); plot(w,angle(H23),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_3) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});

figure(33)
set(gcf,'position',[70, 70, 1400, 500])
subplot(2,3,1); plot(w,abs(H31),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [rad/N]'); title('H_3_,_1'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,2); plot(w,abs(H32),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [rad/N]'); title('H_3_,_2'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,3); plot(w,abs(H33),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [rad/N]'); title('H_3_,_3'); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,4); plot(w,angle(H31),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_1) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,5); plot(w,angle(H32),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_2) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});
subplot(2,3,6); plot(w,angle(H33),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_3) [deg]'); yticks(y_ticks); ylim(y_lim); xline([w_d_ord(1), w_d_ord(3), w_d_ord(5)], '-r', {'w1', 'w2', 'w3'});

%% stessa figura dividendola in 2: amplitudes da phases
figure(41)
set(gcf,'position',[200, 70, 1200, 700])
subplot(3,3,1); plot(w/2/pi,abs(H11),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_,_1');
subplot(3,3,2); plot(w/2/pi,abs(H12),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_,_2');
subplot(3,3,3); plot(w/2/pi,abs(H13),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_,_3');

subplot(3,3,4); plot(w/2/pi,abs(H21),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [rad/N]'); title('H_2_,_1');
subplot(3,3,5); plot(w/2/pi,abs(H22),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [rad/N]'); title('H_2_,_2');
subplot(3,3,6); plot(w/2/pi,abs(H23),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_3| [rad/N]'); title('H_2_,_3');

subplot(3,3,7); plot(w/2/pi,abs(H31),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [rad/N]'); title('H_3_,_1');
subplot(3,3,8); plot(w/2/pi,abs(H32),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [rad/N]'); title('H_3_,_2');
subplot(3,3,9); plot(w/2/pi,abs(H33),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [rad/N]'); title('H_3_,_3');

figure(42)
set(gcf,'position',[200, 70, 1200, 700])
subplot(3,3,1); plot(w/2/pi,angle(H11)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_1) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_1_,_1');
subplot(3,3,2); plot(w/2/pi,angle(H12)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_2) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_1_,_2');
subplot(3,3,3); plot(w/2/pi,angle(H13)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_3) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_1_,_3');

subplot(3,3,4); plot(w/2/pi,angle(H21)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_1) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_2_,_1');
subplot(3,3,5); plot(w/2/pi,angle(H22)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_2) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_2_,_2');
subplot(3,3,6); plot(w/2/pi,angle(H23)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_3) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_2_,_3');

subplot(3,3,7); plot(w/2/pi,angle(H31)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_1) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_3_,_1');
subplot(3,3,8); plot(w/2/pi,angle(H32)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_2) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_3_,_2');
subplot(3,3,9); plot(w/2/pi,angle(H33)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_3) [deg]'); yticks(y_ticks); ylim(y_lim); title('H_3_,_3');
%% 3.b & 3.c Co-located FRF on displacement of point A and rotation of theta3

Phi =[1 0 0
      0 1 0
      1 0 r];

M_1 = (inv(Phi))'*M_gen*inv(Phi);
K_1 = (inv(Phi))'*K_gen*inv(Phi);
C_1 = (inv(Phi))'*C_gen_ray*inv(Phi);

for ii=1:length(w)
    H_1 = inv(-w(ii)^2*M_1+1i*w(ii)*C_1+K_1);
    H11_1(ii)=H_1(1,1);
    H12_1(ii)=H_1(1,2);
    H13_1(ii)=H_1(1,3);
    H21_1(ii)=H_1(2,1);
    H22_1(ii)=H_1(2,2);
    H23_1(ii)=H_1(2,3);
    H31_1(ii)=H_1(3,1);
    H32_1(ii)=H_1(3,2);
    H33_1(ii)=H_1(3,3);

end

figure(4)
% set(gcf,'position',[200, 70, 1200, 700])
% sgtitle('Frequency Responce Function Amplitudes')

subplot(2,1,1); plot(w/2/pi,abs(H33_1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H''_3_3| [rad/N]'); title('Co-located FRF on point A');
subplot(2,1,2); plot(w/2/pi,angle(H33_1)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H''_3_3) [deg]'); yticks(y_ticks);

figure(51)

subplot(2,1,1); plot(w/2/pi,abs(H33),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [rad/N]'); title(('Co-located FRF on rotation \theta_3'));
subplot(2,1,2); plot(w/2/pi,angle(H33)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_3) [deg]'); yticks(y_ticks);

%%
figure(5)
subplot(3,1,1); plot(w/2/pi,abs(H11_1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H''_1_1| [rad/N]'); title('H''_1_,_1');
subplot(3,1,2); plot(w/2/pi,abs(H22_1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H''_2_2| [rad/N]'); title('H''_2_,_2');
subplot(3,1,3); plot(w/2/pi,abs(H33_1),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H''_3_3| [rad/N]'); title('H''_3_,_3');

%% 3.d Complete time response
t=0:0.01:40;
F = 15*cos(2*pi*0.5*t)+7*cos(2*pi*1.25*t);

jac_F=[1; 0; r];

% H_F_1 = inv(-0.5^2*M_gen+1i*0.5*C_gen_ray+K_gen);
% H_F_2 = inv(-1.25^2*M_gen+1i*1.25*C_gen_ray+K_gen);

H_F_1 = inv((-0.5*2*pi)^2*M_gen+1i*(0.5*2*pi)*C_gen_ray+K_gen);
H_F_2 = inv((-1.25*2*pi)*M_gen+1i*(1.25*2*pi)*C_gen_ray+K_gen);

sum_H_1 = 15*H_F_1*jac_F + 7*H_F_2*jac_F;
sum_H_2 = 0.5*15*H_F_1*jac_F + 1.25*7*H_F_2*jac_F;

const_com = [V_damp_norm_1;lambda_1_ord.'.*V_damp_norm_1]\[pos_t0 - sum_H_1 ;vel_t0 - 1i*sum_H_2];

for ii = 1:length(t)
    free_motion_complete(:,ii) = (const_com.'.*V_damp_norm_1)*exp(lambda_1_ord*t(ii));
end

for ii = 1:length(t)
    steady_state(:,ii) = 15*H_F_1*jac_F*exp(1i*0.5*t(ii)) + 7*H_F_2*jac_F*exp(1i*1.25*t(ii));
end

complete_res = free_motion_complete+steady_state;
figure(7)
set(gcf,'position',[200, 70, 1200, 700])

subplot(4,1,1); plot(t(1,1:1001),F(1,1:1001)); grid on; xlabel('time'); ylabel('amplitude [N]'); title('Forcing term');
subplot(4,1,2); plot(t(1,1:2001), free_motion_complete(:,1:2001), LineWidth=1.2); grid on;legend('y_1(t)','\theta_2(t)', '\theta_3(t)'); xlabel('time');  ylabel({'displacement [m]','rotation [rad]'}); title('Free time response');
subplot(4,1,3); plot(t, steady_state, LineWidth=1.2); grid on; legend('y_1(t)','\theta_2(t)', '\theta_3(t)'); xlabel('time');  ylabel({'displacement [m]','rotation [rad]'}); title('Steady-state response');
subplot(4,1,4); plot(t(1,1:2001), real(complete_res(:,1:2001)), LineWidth=1.2); grid on; ylabel({'displacement [m]','rotation [rad]'}); xlabel('time');  legend('y_1(t)','\theta_2(t)', '\theta_3(t)'); title('Complete time response')


%% 3.e
t=0:0.01:20;

Fdt = 0;
Xdt = 0;

f0 = 0.45;

for ii = 1:5 %number of harmonics used to reconstruct the square wave
    n = 2*ii-1;
    F(ii) = 8 / (pi^2) * (-1)^(ii-1) / n^2;
    w_f(ii) = 2*pi*f0*n; %omega force
    X(:,ii) = inv(-(w_f(ii))^2*M_gen + 1i*w_f(ii)*C_gen_ray + K_gen)*jac_F*F(ii);
    
    Fdt = Fdt + F(ii)*sin(w_f(ii)*t); %force (multiplication for numerator)
    Xdt = Xdt+abs(X(:,ii)).*sin(w_f(ii)*t+angle(X(:,ii)));
end

X_ss_tr =jac_F'*Xdt; %steady state response with triangular force

figure(7)

subplot(2,1,2);
plot(t,X_ss_tr,'linewidth',2); grid on ;xlabel('Time [s]');ylabel('Displacement [m]')
title('Response of y_A')

subplot(2,1,1)
plot(t,Fdt,'k','linewidth',1.5);grid on ;xlabel('Time [s]');ylabel('Amplitude [N]')
title('Forcing term')

%% 4.a

Phi_q = [V_damp_norm_1(:,1), V_damp_norm_1(:,3), V_damp_norm_1(:,5)];
M_q = Phi_q'*M_gen*Phi_q;
K_q = Phi_q'*K_gen*Phi_q;
C_q = Phi_q'*C_gen_ray*Phi_q;

for ii = 1:3
    for jj = 1:3
        if ii~=jj
            M_q(ii,jj)=0;
            K_q(ii,jj)=0;
            C_q(ii,jj)=0;
        end
    end
end

for ii=1:length(w)
    H_mod = inv(-w(ii)^2*M_q+1i*w(ii)*C_q+K_q);
    H11_mod(ii)=H_mod(1,1); H12_mod(ii)=H_mod(1,2); H13_mod(ii)=H_mod(1,3);
    H21_mod(ii)=H_mod(2,1); H22_mod(ii)=H_mod(2,2); H23_mod(ii)=H_mod(2,3);
    H31_mod(ii)=H_mod(3,1); H32_mod(ii)=H_mod(3,2); H33_mod(ii)=H_mod(3,3);
end

figure(8)
set(gcf,'position',[200, 70, 1200, 700])
sgtitle('Frequency Responce Function Amplitudes')
subplot(6,3,1); plot(w/2/pi,abs(H11_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_1| [m/N]'); title('H_1_,_1');
subplot(6,3,2); plot(w/2/pi,abs(H12_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_2| [m/N]'); title('H_1_,_2');
subplot(6,3,3); plot(w/2/pi,abs(H13_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_1_3| [m/N]'); title('H_1_,_3');

subplot(6,3,4); plot(w/2/pi,angle(H11_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_1) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,5); plot(w/2/pi,angle(H12_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_2) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,6); plot(w/2/pi,angle(H13_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_1_3) [rad]'); yticks(y_ticks); ylim(y_lim);

subplot(6,3,7); plot(w/2/pi,abs(H21_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_1| [rad/N]'); title('H_2_,_1');
subplot(6,3,8); plot(w/2/pi,abs(H22_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_2| [rad/N]'); title('H_2_,_2');
subplot(6,3,9); plot(w/2/pi,abs(H23_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_2_13| [rad/N]'); title('H_2_,_3');

subplot(6,3,10); plot(w/2/pi,angle(H21_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_1) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,11); plot(w/2/pi,angle(H22_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_2) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,12); plot(w/2/pi,angle(H23_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_2_3) [rad]'); yticks(y_ticks); ylim(y_lim);

subplot(6,3,13); plot(w/2/pi,abs(H31_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_1| [rad/N]'); title('H_3_,_1');
subplot(6,3,14); plot(w/2/pi,abs(H32_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_2| [rad/N]'); title('H_3_,_2');
subplot(6,3,15); plot(w/2/pi,abs(H33_mod),'b'); grid on; xlabel('Frequency [Hz]'); ylabel('|H_3_3| [rad/N]'); title('H_3_,_3');

subplot(6,3,16); plot(w/2/pi,angle(H31_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_1) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,17); plot(w/2/pi,angle(H32_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_2) [rad]'); yticks(y_ticks); ylim(y_lim);
subplot(6,3,18); plot(w/2/pi,angle(H33_mod)*180/pi,'b'); grid on; xlabel('Frequency [Hz]'); ylabel('\angle(H_3_3) [rad]'); yticks(y_ticks); ylim(y_lim);

%% 4.b co_located FRF of vertical displacement in A with modal approach

FRF_q = [H11_mod;
        H22_mod;
        H33_mod;];


H11_col = 0;
H12_col = 0;
H13_col = 0;
H21_col = 0;
H22_col = 0;
H23_col = 0;
H31_col = 0;
H32_col = 0;
H33_col = 0;
for ii = 1:3 
        H = FRF_q(ii, :);
        H11_col = H11_col + Phi_q(1,ii)*Phi_q(1,ii)*H;
        H12_col = H12_col + Phi_q(1,ii)*Phi_q(2,ii)*H;
        H13_col = H13_col + Phi_q(1,ii)*Phi_q(3,ii)*H;
        H21_col = H21_col + Phi_q(2,ii)*Phi_q(1,ii)*H;
        H22_col = H22_col + Phi_q(2,ii)*Phi_q(2,ii)*H;
        H23_col = H23_col + Phi_q(2,ii)*Phi_q(3,ii)*H;
        H31_col = H31_col + Phi_q(3,ii)*Phi_q(1,ii)*H;
        H32_col = H32_col + Phi_q(3,ii)*Phi_q(2,ii)*H;
        H33_col = H33_col + Phi_q(3,ii)*Phi_q(3,ii)*H;
        H33_col_store(ii,:) = Phi_q(3,ii)*Phi_q(3,ii)*H;
        H33_col_store_A(ii,:) =  Phi_q(1,ii)*Phi_q(1,ii)*H + r*Phi_q(1,ii)*Phi_q(3,ii)*H + r*Phi_q(3,ii)*Phi_q(1,ii)*H + (r^2)*Phi_q(3,ii)*Phi_q(3,ii)*H;
end

H_col_A = H11_col + r*H13_col + r*H31_col + r^2*H33_col;


figure(9)
subplot(2,1,1);
plot(w/2/pi,abs(H33_col_store_A(1,:)),w/2/pi,abs(H33_col_store_A(2,:)),w/2/pi,abs(H33_col_store_A(3,:)),'LineWidth',1.5);
grid on;
hold on;
plot(w/2/pi,abs(H33_1),'--k','LineWidth',3); 
plot(w/2/pi,abs(H_col_A),'LineWidth',1.5);
xlabel('Frequency [Hz]'); ylabel('Amplitude [m/N]')
legend('first mode contribution','second mode contribution','third mode contribution','Co-located |FRF_A| from original coordinates','Co-located |FRF_A| from modal coordinates');

subplot(2,1,2);
plot(w/2/pi,angle(H33_col_store_A(1, :))*180/pi,w/2/pi,angle(H33_col_store_A(2, :))*180/pi,w/2/pi,angle(H33_col_store_A(3, :))*180/pi,'LineWidth',1.5);
grid on;
hold on
plot(w/2/pi,angle(H33_1)*180/pi,'--k','LineWidth',3); 
plot(w/2/pi,angle(H_col_A)*180/pi,'LineWidth',1.5); 
xlabel('Frequency [Hz]');ylabel('Phase [deg]');
legend('first mode contribution','second mode contribution','third mode contribution','Co-located \angle(FRF_A) from original coordinates','Co-located \angle(FRF_A) from modal coordinates');
sgtitle('co-located FRF on point A');
%%
figure(60)
plot(w/2/pi,abs(H33_col_store_A(1,:)))
%% 4.c co_located FRF of rotation theta3 with modal approach

figure(10)
subplot(2,1,1);

plot(w/2/pi,abs(H33_col_store(1,:)),w/2/pi,abs(H33_col_store(2,:)),w/2/pi,abs(H33_col_store(3,:)),'LineWidth',1.5);
grid on;
hold on;
plot(w/2/pi,abs(H33),'--k','LineWidth',3); 
plot(w/2/pi,abs(H33_col),'LineWidth',1.5); 
xlabel('Frequency [Hz]');ylabel('Amplitude [m/N]')
legend('first mode contribution','second mode contribution','third mode contribution','Co-located |FRF_\theta_3| from original coordinates','Co-located |FRF_\theta_3| from modal coordinates');

subplot(2,1,2);

plot(w/2/pi,angle(H33_col_store(1,:))*180/pi,w/2/pi,angle(H33_col_store(2,:))*180/pi,w/2/pi,angle(H33_col_store(3,:))*180/pi,'LineWidth',1.5);
grid on;
hold on;
plot(w/2/pi,angle(H33)*180/pi,'--k','LineWidth',3); 
plot(w/2/pi,angle(H33_col)*180/pi,'LineWidth',1.5);
xlabel('Frequency [Hz]');ylabel('Phase [deg]');
legend('first mode contribution','second mode contribution','third mode contribution','Co-located \angle(FRF_\theta_3) from original coordinates','Co-located \angle(FRF_\theta_3) from modal coordinates');
sgtitle('co-located FRF of rotation \theta_3');

%% 4.d commplete time 

lambda_A_q = Phi_q'*jac_F;

H_A1 = FRF_q(1, :) * lambda_A_q(1);     
H_A2 = FRF_q(2, :) * lambda_A_q(2);
H_A3 = FRF_q(3, :) * lambda_A_q(3);

% single mode (first one)
Phi_q_sm = [Phi_q(:,1), zeros(3,1), zeros(3,1)];

lambda_A_q_sm = Phi_q_sm'*jac_F;

H_A1_sm = FRF_q(1, :) * lambda_A_q_sm(1);     
H_A2_sm = FRF_q(2, :) * lambda_A_q_sm(2);
H_A3_sm = FRF_q(3, :) * lambda_A_q_sm(3);

% Force with f1=0.5 & A1=15
f1=0.5;
A1=15;

omega_force = 2*pi*f1;
index = find(abs(w-omega_force)==min(abs(w-omega_force)));

q1 = abs(H_A1(index)) * A1 * cos(2*pi*f1*t + angle(H_A1(index)));
q2 = abs(H_A2(index)) * A1 * cos(2*pi*f1*t + angle(H_A2(index)));
q3 = abs(H_A3(index)) * A1 * cos(2*pi*f1*t + angle(H_A3(index)));

q1_sm = abs(H_A1_sm(index)) * A1 * cos(2*pi*f1*t + angle(H_A1_sm(index)));
q2_sm = abs(H_A2_sm(index)) * A1 * cos(2*pi*f1*t + angle(H_A2_sm(index)));
q3_sm = abs(H_A3_sm(index)) * A1 * cos(2*pi*f1*t + angle(H_A3_sm(index)));

q = [q1; q2; q3];
q_sm = [q1_sm;q2_sm;q3_sm];

x = Phi_q * q;
x_sm =Phi_q * q_sm;

figure(11)
sgtitle('Steady-state response to F(t) = A_1cos(2\pif_1t)')
subplot (4,1,1);
plot(t,q,'LineWidth',1.2);
title('Complete system response in Modal coordinates');
xlabel('Time [s]');ylabel('displacement');
grid on
legend('q_1(t)','q_2(t)','q_3(t)')

subplot (4,1,2);
plot(t,q_sm,'LineWidth',1.2);
title('First mode response in Modal coordinates');
xlabel('Time [s]');ylabel('displacement');
grid on
legend('q_1(t)','q_2(t)','q_3(t)');

subplot (4,1,3);
plot(t,x,'LineWidth',1.2);
title('Complete system response in Original coordinates');
xlabel('Time [s]');ylabel({'displacement [m]','rotation [rad]'});
grid on
legend('y_1(t)','\theta_2(t)','\theta_3(t)')

subplot (4,1,4);
plot(t,x_sm,'LineWidth',1.2);
title('First mode response in Original coordinates');
xlabel('Time [s]');ylabel({'displacement [m]','rotation [rad]'});
grid on
legend('y_1(t)','\theta_2(t)','\theta_3(t)')




% Force with f2=1.25 & A2=7
f2=1.25;
A1=7;

omega_force_2 = 2*pi*f2;
index_2 = find(abs(w-omega_force_2)==min(abs(w-omega_force_2)));

q1_2 = abs(H_A1(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A1(index_2)));
q2_2 = abs(H_A2(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A2(index_2)));
q3_2 = abs(H_A3(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A3(index_2)));

q1_sm_2 = abs(H_A1_sm(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A1_sm(index_2)));
q2_sm_2 = abs(H_A2_sm(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A2_sm(index_2)));
q3_sm_2 = abs(H_A3_sm(index_2)) * A1 * cos(2*pi*f2*t + angle(H_A3_sm(index_2)));

q_2 = [q1_2; q2_2; q3_2];
q_sm_2 = [q1_sm_2;q2_sm_2;q3_sm_2];

x_2 = Phi_q * q_2;
x_sm_2 =Phi_q * q_sm_2;

figure(12)
sgtitle('Steady-state response to F(t) = A_2cos(2\pif_2t)')
subplot (4,1,1);
plot(t,q_2,'LineWidth',1.2);
title('Complete system response in Modal coordinates');
xlabel('Time [s]');ylabel('displacement');
grid on
legend('q_1(t)','q_2(t)','q_3(t)')

subplot (4,1,2);
plot(t,q_sm_2,'LineWidth',1.2);
title('First mode response in Modal coordinates');
xlabel('Time [s]');ylabel('displacement');
grid on
legend('q_1(t)','q_2(t)','q_3(t)');

subplot (4,1,3);
plot(t,x_2,'LineWidth',1.2);
title('Complete system response in Original coordinates');
xlabel('Time [s]');ylabel({'displacement [m]','rotation [rad]'});
grid on
legend('y_1(t)','\theta_2(t)','\theta_3(t)')

subplot (4,1,4);
plot(t,x_sm_2,'LineWidth',1.2);
title('First mode response in Original coordinates');
xlabel('Time [s]');ylabel({'displacement [m]','rotation [rad]'});
grid on
legend('y_1(t)','\theta_2(t)','\theta_3(t)')

