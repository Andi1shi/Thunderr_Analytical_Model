function [vxt,vyt,vt,alphat] = THUNDERR_PLOT(X,BestFitness,ti,dt,tf,xp,yp,vrt,alphart,fs,xc0_lb,xc0_ub,yc0_lb,yc0_ub,FontSize)

% THUNDERR_PLOT: Analytical Model for Wind Velocity Reconstruction
%
% This function provides a detailed reconstruction of the horizontal mean wind velocity 
% during a downburst event, based on the Teaching Learning Based Optimization (TLBO)
% algorithm's outcomes. It visualizes the effect of a downburst (thunderstorm) on wind 
% velocity using given parameters.
%
% The function also demonstrates how the NMSE (Normalized Mean Square Error) can be used 
% to evaluate the accuracy of the downburst simulation against actual recorded data.
%
% Input Arguments:
% X        - Optimized parameter vector from TLBO algorithm.
% BestFitness - Best fitness value from the optimization, used for analysis.
% ti, tf   - Initial and final time for the simulation.
% dt       - Time step for the simulation.
% xp, yp   - Coordinates representing the position (x, y) of the anemometer in the simulation space.
% vrt, alphart - Recorded wind velocities and directions, used for comparison.
% fs       - Sampling frequency of the recorded data.
% xc0_lb, xc0_ub, yc0_lb, yc0_ub - Bounds for the downburst touch down location.
% FontSize - Font size for plotting the results.
%
% Output Arguments:
% vxt, vyt - Horizontal components of the mean wind velocity.
% vt       - Total mean wind velocity.
% alphat   - Direction of the wind velocity.
%
% Model External Parameters (Derived from X):
% xc0  - x-component of downburst touch down (m).
% yc0  - y-component of downburst touch down (m).
% R    - Downburst downdraft radius (m).
% v1   - Maximum radial velocity in the downburst (m/s).
% T1   - Period of linear intensification (s).
% Tf   - Total duration of the downburst (s).
% v2   - Downburst translation velocity (m/s).
% beta2- Downburst translation direction (degrees from the East).
% v3   - ABL background wind speed (m/s).
% beta3- ABL background wind direction (degrees from the East).

xc0 = X(1);        % - x  component of  downburst touch down (m)
yc0 = X(2);        % - y  component of  downburst touch down (m)
R = X(3);          % - downburst downdraft radius (m)
psi = X(4);        % - psi = Rmax/R Radius of maximum radial velocity 
v1 =X(5);          % - downburst maximum radial velocity (m/s)
T1 = X(6);         % - period of linear intensification of the downburst (s)
Tf = X(7);         % - total duration of the downburst (s)
v2 = X(8);         % - downburst translation velocity (m/s)
beta2 = X(9);      % - downburst translation direction (polar angle from the East) (deg)
v3 = X(10);        % - ABL background wind speed (m/s)
beta3 = X(11);     % - ABL background wind direction (polar angle from the East) (deg)
%-------------------------------------------------------------------------%

%% -------------------Model Necessary Parameters -------------------------%
% ti       - initial simulation time (s), usually = 0 s
% dt       - simulation time step (s)
% tf       - final simulation time (s)
% xp       - x station location component 
% yp       - y station location component
% vxrt     - moving average(T_ma = 30 s) recorded (r) data along x - axis
% vyrt     - moving average(T_ma = 30 s) recorded (r) data along y - axis
% vrt      - moving average(T_ma = 30 s) wind speed recorded (r) data 
% alphart  - moving average(T_ma = 30 s) wind direction recorded (r) data
% fs       - recorded data sampling frequency (10 Hz mostly)
%-------------------------------------------------------------------------%
%% -------------------Model Internal Parameters----------------------------
% Internal parameters -  Holmes and Oliver function 2000
Rs = R;                 % Exponential decay constant (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rmax = 2*R;             % Radius of maximum radial velocity (m)
Rmax = psi*R;             % Radius of maximum radial velocity (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Internal parameters -  Auxiliary function translational component
RE = 2*Rmax;             % radius of maximum outflow extension (m)
a = 0.2;                 % smooth factor
RE_max = (1 + a)*RE;
%-------------------------------------------------------------------------%
% Internal parameters -  Intensity decay function translational component
b = 0.2;                 % smooth factor     
Ti = b.*T1;                    
Tf_max = (1 + b).*Tf;  
%-------------------------------------------------------------------------%
%% Starting of Thunderr_02
% Time axis
t = (ti:dt:tf-dt)';
nt = dt*fs;

%% ----------------------(1) - Stationary downburst-----------------------%

% Holmes and Oliver function 2000

vr = @(s) (v1.*s./Rmax).*(s <= Rmax)...
                      + (v1.*exp(-(((s-Rmax)./Rs).^2))).*(s > Rmax);
% Note
% s is a dummy variable
%-------------------------------------------------------------------------%

% Chay at al. 2006 - intensity decat function radial component
% Enssuring that Tf is grater than T1
if Tf <= T1           
    error (' Tf <= T1, Tf shoud be > T1')
end
% Exponential decay time constant (sec)
c = (Tf - T1)./(log(10));         
PIr = @(s) (s./T1).*(s<=T1) + (exp(-((s-T1)./c))).*(s>T1);
%-------------------------------------------------------------------------%

%% ----------------------(2) - Translating downburst----------------------%

% Downburst path
xct = xc0 + v2*cosd(beta2)*t;
yct = yc0 + v2*sind(beta2)*t;
%-------------------------------------------------------------------------%
% Relative distance between the downburst center and the station
rt = sqrt(((xp - xct).^2 + (yp - yct).^2));

% Relative direction between the downburst center and the station
% (polar angle i.e. from East)
beta1t = atan2d((yp - yct),(xp - xct));
%-------------------------------------------------------------------------%
% Xhelaj et al. (2020) Auxiliary function translational component

Delta = @(s) (1).*((0 <= s)&(s <= RE))...
              + ((1/2).*(1 + cos((pi/(a*RE)).*(s - RE)))).*((RE < s)&(s <= RE_max))...
              + (0).*(s > RE_max);
%-------------------------------------------------------------------------%
% Xhelaj et al. (2020) Intensity - decay function translational component

PIt =  @(s) ((1/2).*(1 + cos((pi./Ti).*(s - Ti)))).*((0 <= s)&(s <= Ti))...
          + (1).*((Ti < s)&(s <= Tf))...
          + ((1/2).*(1 + cos((pi./(b.*Tf)).*(s - Tf)))).*((Tf < s)&(s <= Tf_max))...
          + (0).*(s > Tf_max);

%% ------------------------ Simulated data---------------------------------
vxt = vr(rt).*PIr(t).*cosd(beta1t) + v2.*Delta(rt).*PIt(t).*cosd(beta2) + v3*cosd(beta3);
vyt = vr(rt).*PIr(t).*sind(beta1t) + v2.*Delta(rt).*PIt(t).*sind(beta2) + v3*sind(beta3);

vt = sqrt(vxt.^2 + vyt.^2);
betat = atan2d(vyt,vxt);
alphat = mod(270 - betat,360);
%--------------------------------------------------------------------------
%% -------------------------- Recorded data--------------------------------

vrt = vrt(1:nt:end);
alphart = alphart(1:nt:end);
%--------------------------------------------------------------------------
%% --------------------- Shifting the  simulated data ---------------------
% Shifting the  data
[~,ind1] = max(vrt);
[~,ind2] = max(vt);
Shift_abs = abs(ind2 - ind1);
if  ind1 > ind2
    Shift = Shift_abs;
else
    Shift = - Shift_abs;
end


vt = circshift(vt,Shift);
alphat = circshift(alphat,Shift);

%% -----------------------------Plotting ----------------------------------
% Downdraft radius
ksi = 0:0.01:2*pi;    % dummy variable
xd = xc0 + R*cos(ksi);
yd = yc0 + R*sin(ksi);

% Max outflow velocit
xv_max = xc0 + Rmax*cos(ksi);
yv_max = yc0 + Rmax*sin(ksi);
% Outflow extension
xout = xc0 + RE_max*cos(ksi);
yout = yc0 + RE_max*sin(ksi);

fig = figure;
fig.WindowState = 'maximized';
subplot(2,2,[1 3])
plot(xct,yct,'--b','LineWidth',2);
hold on
p0 = plot(xc0,yc0,'ok','MarkerSize',8);
plot([xc0_lb xc0_ub],[0 0],'--k','LineWidth',2);
plot([0 0],[yc0_lb yc0_ub],'--k','LineWidth',2)
p1 = plot(xd,yd,'r','LineWidth',2);
p2 = plot(xv_max,yv_max,'--r','LineWidth',2);
p3 = plot(xout,yout,'b','LineWidth',2);
p4 = plot(xp,yp,'ok','MarkerFaceColor','k','MarkerSize',8);


str1 = ['x_{c0} = ',num2str(round(xc0,2)),' m'];
str2 = ['y_{c0} = ',num2str(round(yc0,2)),' m'];
str3 = ['R = ',num2str(round(R,2)),' m'];
str4 = ['R_{max} = ',num2str(round(Rmax,2)),' m'];
str5 = ['v_{r_{max}} = ',num2str(round(v1,2)),' m/s'];
T1 = T1/60;
str6 = ['T_1 = ',num2str(round(T1,2)),' min'];
Tf = Tf/60;
str7 = ['T_f = ',num2str(round(Tf,2)),' min'];
str8 = ['v_{trans} = ',num2str(round(v2,2)),' m/s'];
alpha2 = mod(270 - beta2,360);
str9 = ['alpha_{trans} = ',num2str(round(alpha2,2)),' deg'];
str10 = ['v_{back} = ',num2str(round(v3,2)),' m/s'];
alpha3 = mod(270 - beta3,360);
str11 = ['alpha_{back} = ',num2str(round(alpha3,2)),' deg'];

tx = xc0_ub;
ty = yc0_ub;
dt =  1000;
text(tx, ty,str1,'Color','k','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-dt,str2,'Color','k','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-2*dt,str3,'Color','r','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-3*dt,str4,'Color','r','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-4*dt,str5,'Color','r','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-5*dt,str6,'Color','b','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-6*dt,str7,'Color','b','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-7*dt,str8,'Color','b','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-8*dt,str9,'Color','b','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-9*dt,str10,'Color','m','FontSize',12,'HorizontalAlignment','left')
text(tx, ty-10*dt,str11,'Color','m','FontSize',12,'HorizontalAlignment','left')

xlabel('x (m)')
ylabel('y (m)')
title(['ALGORITHM: STLBO     Best Fitness Value = ',num2str(BestFitness)],'Color','k')
grid minor
legend([p0 p1 p2 p3 p4],'touch down','downburst downdraft','downburst max outflow','downburst transalation region','anemometer')
axis equal
xlim([xc0_lb xc0_ub]);
ylim([yc0_lb yc0_ub]);
ax = gca; ax. FontSize = FontSize;

subplot(2,2,2)
plot(t,vt,'.-r')
xlabel('t (s)')
ylabel('v (m/s)')
hold on
plot(t,vrt,'.-b')
xlim([t(1) t(end)])
ylim([0 max(vrt)+2])
grid on
ax = gca; ax. FontSize = FontSize;

subplot(2,2,4)
plot(t,alphat,'.r')
hold on
plot(t,alphart,'.b')
xlim([t(1) t(end)])
ylim([0 360])
yticks(0:45:360)
xlabel('t (s)')
ylabel('\alpha (t) (deg)')
grid on
ax = gca; ax. FontSize = FontSize;
%-------------------------------------------------------------------------%

end