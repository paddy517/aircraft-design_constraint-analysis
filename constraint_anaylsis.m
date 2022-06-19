%% cruising equation

% parameters collection
W_to = 260000 ; % take off weight(kg)
W = 0.78*W_to ; % current weight
p = 0.3428 ;% (kg/m^3)
afa = 0.5 ;% engine efficiency
V_cruising = 257.25; %crusing velocity (m/s)
k1 = 0.006027; 
cd_0 = 0.0085;
q = 0.5*p*V_cruising^2 ;% dynamic presuure
s = 370;% wing surface (m^2)
beta = W/W_to;
W_to_s = linspace(10,800,400);
T_W_cruise = (beta/afa)*(k1*(beta/q)*(W_to_s)+(cd_0*q /beta)*(1./W_to_s));  % thrust loading for crusing

plot(W_to_s,T_W_cruise)
hold on

%% climb  equation

% parameters collection
W = 0.96*W_to ;% weight
afa = 0.5; % engine efficiency
k1 = 0.006027;
cd_0 = 0.0085;
V_climb = 154.3 ; % velocity of climb m/s
q = 0.5*p*V_climb^2 ; % dynamic presuure
s = 370 ;% wing surface (m^2)
dh_dt = 10 ;% climb rate m/s
beta = W/W_to;
W_to_s = linspace(10,800,400);
T_W_climb = (beta/afa)*(k1*(beta/q)*(W_to_s)+(cd_0*q/(beta)*(1./W_to_s))+(1/V_climb)*(dh_dt)); % thrust loading for climb aand decent

plot(W_to_s,T_W_climb)
hold on

%% loister / speed turn
W = 0.61*W_to;% weight
afa = 0.6175; % engine efficiency
k1 = 0.006027;
cd_0 = 0.0085;
V_loister = 154; %velocity of loister m/s
q = 0.5*p*V_loister^2; % dynamic presuure
s = 370 ;% wing surface m^2
V_lositer_knots = 300 ;% velocity knots
zeta = (V_lositer_knots/10) + 7 ; % unit V: knots 
n = 1/sind(zeta);
beta = W/W_to; 
W_to_s = linspace(10,800,400);
T_W_loister = (beta/afa)*(k1*n^2*(beta/q)*(W_to_s)+(cd_0*q/beta)*(1./W_to_s));  % thrust loading for loister

plot(W_to_s,T_W_loister)
hold on

%% take off ground roll
W =  260000; % weight kg
p =  1.225; % air density  kg/m3
V_stall = 82; % m/s
Cl_max = 2.8;
s =  370; % wing surface m2
beta = 0.99;
W_to_s = linspace(10,800,400);
T_W_takeoff = sqrt((2*beta)*(W_to_s)/(p*(V_stall)^2*Cl_max)); % thrust loading for take off 

plot(W_to_s,T_W_takeoff)
hold on


%% stall speed
W = 260000; % weight (kg)
s = 370; % wing surface
g= 9.8;
Cl_max = 2.8;
p = 1.225;% air density
V_stall = 82.42; % m/s %((2*W*g)/(p*s*Cl_max))^0.5 %1.3*vs;
T_Wmax = 30;
W_to_s = 0.5*p*V_stall^2*Cl_max /(g*1.3^2); % weight/surface

plot([W_to_s,W_to_s],[-0.01,T_Wmax])
hold on

%% plot

legend('cruising','climb','loister','take off ground roll','stall speed')
xlabel('W/S (kg/m^2)');
ylabel('T/W ');