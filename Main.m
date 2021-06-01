clc; clear; close all;
Input_menu; % Main calls Input_menu
load('_data.mat');

%% Data input and treatment

% Ambient
SIM(1) = data.ambient(1); %K (Tamb)
SIM(2) = data.ambient(2)*10^5; %bar > Pa (Pamb)

% Simulator
SIM(3) = data.simulation(1); %s (dt_sim)
SIM(4) = data.simulation(2); %(-) (N_sim)
SIM(5) = data.simulation(3); %(-) (n_sim) / UNUSED
SIM(6) = data.simulation(4); %s (t_sim)
 
% Chamber - Geometry
CH(1) = data.chamber.geometry(1)*10^(-3); %mm > m (Dext_ch)
CH(2) = data.chamber.geometry(2)*10^(-3); %mm > m (Dint_ch)
CH(3) = data.chamber.geometry(3)*10^(-3); %mm > m (L_ch)
CH(4) = data.chamber.geometry(4); %N (Blk_union)

% Liner & Insulator
CH(5) = data.liner(1)*10^(-3); %mm > m (t_li)
%c_li = data.liner(2); %J/kgK
%lambda_li = data.liner(3); %W/mK

CH(6) = data.insulator(1)*10^(-3); %mm > m (t_ins)
%c_ins = data.insulator(2); %J/kgK
%lambda_ins = data.insulator(3); %W/mK

% Chamber - Material
CH(7) = data.chamber.material.properties(1)*10^9; %GPa > Pa (E_ch)
CH(8) = data.chamber.material.properties(2); %kg/m3 (dens_ch)
CH(9) = data.chamber.material.properties(3); %(-) (poisson_ch)
CH(10) = data.chamber.material.properties(4)+273; %ºC > K (Tmelt_ch)
CH(11) = data.chamber.material.properties(5); %W/mK (lambda_ch)

yield_ch_T = data.chamber.material.yield(1,:)+273; %ºC > K
yield_ch = data.chamber.material.yield(2,:)*10^6; %MPa > Pa

uts_ch_T = data.chamber.material.uts(1,:)+273; %ºC > K
uts_ch = data.chamber.material.uts(2,:)*10^6; %MPa > Pa

% Nozzle - Geometry
NOZ(1) = data.nozzle.geometry(1)*10^(-3); %mm > m (Di_no)
NOZ(2) = data.nozzle.geometry(2)*10^(-3); %mm > m (Dt_no)
NOZ(3) = data.nozzle.geometry(3)*10^(-3); %mm > m (De_no)
NOZ(4) = data.nozzle.geometry(4); %º (conv_no)
NOZ(5) = data.nozzle.geometry(5); %º (div_no)
NOZ(6) = data.nozzle.geometry(6); %N (nozzle_union)

% Nozzle - Material
NOZ(7) = data.nozzle.material.properties(1)*10^9; %GPa > Pa (E_no)
NOZ(8) = data.nozzle.material.properties(2); %kg/m3 (dens_no)
NOZ(9) = data.nozzle.material.properties(3); %(-) (poisson_no) / UNUSED
NOZ(10) = data.nozzle.material.properties(4)+273; %ºC > K (Tmelt_no)
NOZ(11) = data.nozzle.material.properties(5); %W/mK (lambda_no)

yield_no_T = data.nozzle.material.yield(1,:)+273; %ºC > K
yield_no = data.nozzle.material.yield(2,:)*10^6; %MPa > Pa

uts_no_T = data.nozzle.material.uts(1,:)+273; %ºC > K
uts_no = data.nozzle.material.uts(2,:)*10^6; %MPa > Pa

% Propellant - Geometry
PROP(1) = data.propellant.geometry(1)*10^(-3); %mm > m (Dext_p)
PROP(2) = data.propellant.geometry(2)*10^(-3); %mm > m (Dint_p)
PROP(3) = data.propellant.geometry(3)*10^(-3); %mm > m (L0_p)
PROP(4) = data.propellant.geometry(4); %(-) (N_p)
PROP(5) = data.propellant.geometry(5); %(-) (Outside inhibited)
PROP(6) = data.propellant.geometry(6); %(-) (Ends inhibited)

% Propellant - Properties
PROP(7) = data.propellant.properties(1); %mm/sec / MPa (ap)
PROP(8) = data.propellant.properties(2); %(-) (np)
PROP(9) = data.propellant.properties(3); %kg/m3 (dens_p)
PROP(10) = data.propellant.properties(4); %(-) (ratio_p)
PROP(9) = PROP(9)*PROP(10); % Rewrittes PROP(9) as the real density
PROP(11) = data.propellant.properties(5); %W/mK (lambda_p)
PROP(12) = data.propellant.properties(6); %J/kgK (c_p)
PROP(13) = data.propellant.properties(7)*10^(-3); %g/mol > kg/mol (Mw_g)
PROP(14) = data.propellant.properties(8); %J/kgK (Cv_g)
PROP(15) = data.propellant.properties(9); %(-) (k_g)
PROP(16) = data.propellant.properties(10); %Pa·s (mu_g)
PROP(17) = data.propellant.properties(11); %W/mK (lambda_g)
PROP(18) = data.propellant.properties(12); %K (T_flame_p)
PROP(19) = data.propellant.properties(13); %K (T_autoig_p)
PROP(20) = data.propellant.properties(14); %for SI units (alpha erosive burn)
PROP(21) = data.propellant.properties(15); %(-) (beta erosive burn)
PROP(22) = data.propellant.properties(16); %(-) (eff_comb)
PROP(23) = data.propellant.properties(17); %(Xw) % It's hardcoded although only true for A24

%% Simulation time
t_vec = 0:SIM(3):SIM(6); % The maximum duration evaluated. The simulation should stop before
t_vec = t_vec.';

%% FV Discretization
[FV_type,FV_length,FV_area] = FV_discr (CH(3),PROP(3),PROP(4),SIM(4),PROP(1),PROP(2)); 
% Only works for cylindrical and equally sized BATES
FV_area_old=FV_area; 
FV_area = zeros(2,length(FV_type));
FV_dx = zeros(2,length(FV_type));
FV_area(2,:) = FV_area_old;
clear FV_area_old;
FV_dx(2,:) = (FV_length(:,2)-FV_length(:,1)).';
SIM(4) = length(FV_type);  % If discr has created more FV
%% Ignition
FV_ignited = zeros(2,length(FV_type));

FV_Twall = zeros(2,length(FV_type));
FV_Twall(2,:) = SIM(1);

% Additional parameters definitions
PROP(24) = 0.3; % in percentage (pi_k_p) HARDCODED, in truth it oscillates
PROP(24) = PROP(24)*0.01; % in per one (pi_k_p)
PROP(25) = PROP(24)*(1-PROP(8)); % (sigma_p)
PROP(26) = 8.314/(PROP(13)); % (r_g)

PROP(7) = PROP(7)*exp((SIM(1)-293)*PROP(25)); % Reference temperature is considered to be ambient. HARDCODED.
% Recalculation of ap depending on propellant initial temperature

CH(12) = pi/4*CH(2)^2*CH(3); % (V_ch)
CH(13) = pi/4*(CH(2)^2-PROP(1)^2)*CH(3); % (V_liner_ins)

sigma_rad = 5.6703*10^(-8); % Stefan-Boltzmann constant

%% Losses
% Divergence * Boundary layer * Reactions in nozzle * Two phase
eff_tot = 0.5*(1+cos(NOZ(5)*pi/180))*0.985*0.995*0.95; % HARDCODED

ero_f = 1.06; % 6% expected throat increase due to erosion
dDt_no = (NOZ(2)*ero_f-NOZ(2))/t_vec(end);
% This is the the nozzle erosion in m/s. It should be noted, that this
% consider that the burn lasts the full simulated time. This is incorrect
% but it's impossible to estimate this ratio beforehand otherwise

%% Vector initialization
% Dimensions
Dint_p_t = zeros(2,length(FV_type)); % internal diameter as a function of time
Dt_no_t = zeros(length(t_vec),1); % throat erosion as a function of time

% Mass flows 
mdot_gen_el = zeros(1,length(FV_type));

% Gas parameters
T_g_el = zeros(2,length(FV_type)); % Element gas temperature

% Heat exchange
Qw_el = zeros(1,length(FV_type));

% Results
CF = zeros(length(t_vec),1);
cstar = zeros(length(t_vec),1);
T = zeros(length(t_vec),1);
m_p = zeros(length(t_vec),1);
Isp = zeros(length(t_vec),1);

% Additional flags
f_fail = false;

%% Initial conditions
P_ch = zeros(length(t_vec),1);
T_ch = zeros(length(t_vec),1);
P_ch(1) = SIM(2);
T_ch(1) = SIM(1);

Dint_p_t(2,:) = PROP(2);
FV_Twall(2,:) = SIM(1);
T_g_el(2,:) = SIM(1);
Dt_no_t(1) = NOZ(2);



if 1
    fprintf('Select the amount of cracking presented by the propellant.\n')
    cr_type = input('(1 - None / 2 - Microcracks / 3 - Bridged gaps / 4 - Open gaps): ');
    if cr_type~=1 && cr_type~=2 && cr_type~=3 && cr_type~=4
        cr_type = 1;
    end
else
    cr_type = 1;
end
cr_factor = [1 1.28 4 8];
% Something higher than microcracks will break most chambers    


%% Igniters
fprintf('Choose a csv igniter curve to use.\n')
[file,path] = uigetfile('.\Libraries\Ign\*.csv');
fullpath = cat(2,path,file);
ign_read = readtable(fullpath);

IGN = zeros(3,ign_read{1,1});
IGN(1,2) = ign_read{1,2}; % T_ign
IGN(1,3) = ign_read{1,3}; % k
IGN(1,4) = ign_read{1,4}; % Cv
IGN(2,:) = ign_read{2,1:ign_read{1,1}};
IGN(3,:) = ign_read{3,1:ign_read{1,1}};

ximp = mean((FV_length(1,1)+FV_length(1,2))/2); % default position
ximp = round(ximp,3);
fprintf('Introduce the flame impingement point (m)');
ximp = input(['(max: ',num2str(CH(3)),' m / default: ',num2str(ximp),' m): ']);
if ximp>CH(3) || ximp < 0  %Default
    disp('Wrong value. Default ximp is used')
    ximp = mean((FV_length(1,1)+FV_length(1,2))/2); % average position of the ignited 
end

tic; % Initiates chronometer

%% Pressure ratio calculation
At_Ae = NOZ(2)^2/NOZ(3)^2;
for PR = 1:0.01:1000
    ra = ((PROP(15)+1)/2)^(1/(PROP(15)-1))*(PR)^(-1/PROP(15))*sqrt((PROP(15)+1)/(PROP(15)-1)*(1-(PR)^((1-PROP(15))/PROP(15))));
    if abs(At_Ae-ra)/ra < 0.001
        break;
    end
end
% This function finds the ratio between chamber pressure and outlet
% pressure at the nozzle. This value is considered constant

%% CALC: Transient Regime + Steady-State
f_end = false;
if PROP(5) ~= 1 % It doesn't allow execution if the grain exterior is uninhibited
    fprintf('The outside cannot be left uninhibited for this simulator!');
    PROP(5) = 1;
end
it = 1;
k1 = [1 1 1 1 1 1];
k2 = [0 0 0 0 0];
% These are workarounds for ELEMENT SWITCH

it1 = 1; it2 = 4;
% These are workarounds to log into the console the results

dt_save = SIM(3);
while ~f_end
    % Loads previous time step as initial values    
    T_g_el(1,:) = T_g_el(2,:);
    FV_Twall(1,:) = FV_Twall(2,:);
    FV_area(1,:) = FV_area(2,:);
    FV_dx(1,:) = FV_dx(2,:);
    FV_ignited(1,:) = FV_ignited(2,:);
    Dint_p_t(1,:) = Dint_p_t(2,:);
       
    T_flame = flame_temp(PROP(18),PROP(22));
    
    if t_vec(it) <= IGN(2,end)
        mdot_ign = interp1(IGN(2,:), IGN(3,:), t_vec(it));
    else
        mdot_ign = 0; 
    end 
    
    rdot_nom = cr_factor(cr_type)*(PROP(7)*(P_ch(it)*10^(-6))^(PROP(8))*10^(-3));
    % De Vielle equation. For pressures in MPa yields mm/s so they have to be
    % transformed into m/s
    Ab_tot = sum(FV_area(1,FV_ignited(1,:)== 1)); % Full are currently burning
    mdot_gen = PROP(9)*Ab_tot*rdot_nom; % Generated mass flow
    
    if mdot_gen+mdot_ign >= choked_nozzle(Dt_no_t(it),P_ch(it),T_ch(it),SIM(2),PROP(15),PROP(26),1)
        % Choked nozzle if mass flow is high enough
        mdot_noz = choked_nozzle(Dt_no_t(it),P_ch(it),T_ch(it),SIM(2),PROP(15),PROP(26),1);
    else
      	mdot_noz = mdot_gen+mdot_ign;
    end
    V_p = sum(pi/4*(PROP(1)^2-Dint_p_t(1,:).*Dint_p_t(1,:)).*FV_dx(1,:));
    m_p = V_p*PROP(9); % Total propellant mass
    
    V_free = CH(12) - V_p - CH(13); % Chamber without propellant or liner/insulator
    dT_dt = PROP(26)*T_ch(it)/(P_ch(it)*V_free*PROP(15))*(mdot_gen*PROP(15)*PROP(14)...
        *T_flame+mdot_ign*IGN(1,3)*IGN(1,4)*IGN(1,2)-mdot_noz*PROP(15)*PROP(14)*T_ch(it)-...
        T_ch(it)*(PROP(14)*mdot_gen+IGN(1,4)*mdot_ign-PROP(14)*mdot_noz));
    dP_dt =  PROP(26)*T_ch(it)/V_free*(mdot_gen+mdot_ign-mdot_noz)...
        +P_ch(it)/T_ch(it)*dT_dt-P_ch(it)/V_free*Ab_tot*rdot_nom;
    
    for i=1:length(FV_type)-k2(1) % XXXXX START OF ELEMENTWISE XXXXX
        if i==1
            mdot_in_el = mdot_ign;
        else
            mdot_in_el = mdot_out_el; % Carries over from the last iteration
        end
        
        if FV_type(i) == -1 % Jumps burnout elements
            for j=1:length(FV_type)
                if FV_dx(1,j) >= 0 % Elementwise resumed at the next non burnout element
                    i = j; 
                    break;
                else
                end
            end  
        end
        
        V_free_el = pi/4*Dint_p_t(1,i)^2*FV_dx(1,i);
        dens_g_el = P_ch(it)/(PROP(26)*T_g_el(1,i)); % Gas density applying continuity
        v_g_el = mdot_in_el/(dens_g_el*pi/4*Dint_p_t(1,i)^2);
        
        if FV_ignited(1,i) == 1
            rdot_nom_el = PROP(7)*(P_ch(it)*10^(-6))^(PROP(8))*10^(-3);  % ACTUALLY THE SAME AS RDOT_NOM(it)
            rdot_ero_el = erosive_burning(dens_g_el,v_g_el,Dint_p_t(1,i),PROP(20),PROP(21),rdot_nom_el,PROP(9));  
            
            if rdot_ero_el < 0 || ~isreal(rdot_ero_el)
                rdot_ero_el = 0;
            end
            rdot_el = cr_factor(cr_type)*(rdot_nom_el + rdot_ero_el);
            if ~FV_type(i) % Is end element
                mdot_gen_el(i) = PROP(9)*FV_area(1,i)*rdot_el + (PROP(6)-1)*PROP(9)*pi/4*(PROP(1)^2-Dint_p_t(1,i)^2)*rdot_el;
                % Adds the propellant generated by the non-core parts
            else
                mdot_gen_el(i) = PROP(9)*FV_area(1,i)*rdot_el;
                % Only propellant generated by core
            end  
        else
            mdot_gen_el(i) = 0; 
        end
        
        xel = FV_length(i,2); % Position of element end (used as characteristic length)
        Re_g_el = dens_g_el*v_g_el*xel/PROP(16);
        if (t_vec(it) < IGN(2,end)) && (FV_ignited(1,i) == 0)
            if xel < ximp  % Upstream
                h_c_el = 3*Re_g_el^0.8*PROP(17)/Dint_p_t(1,i)*(0.0154+0.125*xel/ximp)*(Dint_p_t(1,i)^2/Dt_no_t(it)^2)^0.4;
            else % Downstream
                h_c_el = 3*Re_g_el^0.8*PROP(17)/Dint_p_t(1,i)*(0.022+0.11*exp(-0.5*(xel-ximp)/Dint_p_t(1,i)))*(Dint_p_t(1,i)^2/Dt_no_t(it)^2)^0.4;
            end
        elseif (t_vec(it) >= IGN(2,end)) && (FV_ignited(1,i) == 0) % Igniter consumed
            h_c_el = Re_g_el^0.8*PROP(17)/Dint_p_t(1,i)*(0.022+0.11*exp(-0.5*xel/Dint_p_t(1,i)))*(Dint_p_t(1,i)^2/Dt_no_t(it)^2)^0.4;
        else % These equations are not contemplated when propellant is burning
            h_c_el = 0;
        end
        
        if FV_ignited(1,i) == 0
            h_r_el = 0.25*sigma_rad*(T_g_el(1,i)^2+FV_Twall(1,i)^2)*(T_g_el(1,i)+FV_Twall(1,i));  
            qw_el = (h_c_el+h_r_el)*(FV_Twall(1,i)-T_g_el(1,i));
            Qw_el(i) = qw_el*FV_area(1,i);
        end   

        
        if FV_ignited(1,i) == 1
            dV_dt_el = FV_area(1,i)*rdot_el; % Like lost volume, will be zeros if propellant doesnt burn
        else
            dV_dt_el = 0;
        end
        
        mdot_stored_el = V_free_el/(PROP(26)*T_g_el(1,i))*(dP_dt-P_ch(it)/T_g_el(1,i)*dT_dt+P_ch(it)/V_free_el*dV_dt_el);
        % Considers dT/dt to be constant
        
        mdot_out_el = mdot_in_el + mdot_gen_el(i) - mdot_stored_el;
        if mdot_out_el < 0
            mdot_out_el = 0;
        end
        
        dT_dt_el = PROP(26)*T_g_el(1,i)/(P_ch(it)*V_free_el*PROP(14))*...
        (mdot_in_el*PROP(15)*PROP(14)*T_g_el(1,i)+mdot_gen_el(i)*PROP(15)*PROP(14)*T_flame-...
        mdot_out_el*T_g_el(1,i)-T_g_el(1,i)*(PROP(14)*mdot_in_el+...
        PROP(14)*mdot_gen_el(i)-PROP(14)*mdot_out_el))+Qw_el(i)/PROP(14);
        
        dP_dt_el = PROP(26)*T_g_el(1,i)/V_free_el*(mdot_in_el+...
        mdot_gen_el(i)-mdot_out_el)-P_ch(it)/V_free_el*dV_dt_el+...
        P_ch(it)/T_g_el(1,i)*dT_dt_el;
        
        if FV_ignited(1,i) == 0
            % Euler for gas elemental temperature
            T_g_el(2,i) = T_g_el(1,i) + dT_dt_el*SIM(3);
            if T_g_el(2,i) > T_ch(it)
                T_g_el(2,i) = T_ch(it); 
                % As the element is part of the chamber it must never surpass
                % chamber temperature
            end
        else
            T_g_el(2,i) = T_ch(it);
        end
        
        if FV_ignited(1,i) == 0 % Element is not ignited
            diff_p = PROP(11)/(PROP(12)*PROP(9));
            h_tot = h_r_el + h_c_el;
            if (FV_Twall(1,i) - SIM(1)) > 10 % Ensure no infinite differential
                dTwall_dt = 4*diff_p*(h_tot/PROP(11)*(T_g_el(1,i)-...
                FV_Twall(1,i)))^3/(3*(2*h_tot/PROP(11)*(T_g_el(1,i)-...
                FV_Twall(1,i))*(FV_Twall(1,i)-SIM(1))+...
                h_tot/PROP(11)*(FV_Twall(1,i)-SIM(1))^2)); % RK4 ignores it, but no matter
                FV_Twall(2,i) = Twall_4th_order_RK (diff_p,h_c_el,h_r_el,PROP(11),SIM(1),FV_Twall(1,i),T_g_el(2,i),T_g_el(1,i),SIM(3)); % Uses RK4
                % Calls 4th Order Runge-Kutta
            else
                dTwall_dt = -qw_el/(FV_area(1,i)*PROP(12)*PROP(9)*0.002*10^(-3));
                FV_Twall(2,i) = FV_Twall(1,i) + dTwall_dt*SIM(3);
                % Euler method for the wall temperature
            end            
            if FV_Twall(2,i) > PROP(19)
                FV_ignited(2,i) = 1;
                disp(['Element ',num2str(i),' is ignited at t = ',num2str(round(t_vec(it),5)),' s']);
                % Ignition message
                FV_Twall(2,i) = T_flame;
            else
                FV_ignited(2,i) = 0;
            end            
        else % Element is ignited
            FV_Twall(2,i) = T_flame;
            FV_ignited(2,i) = 1; % Reafirms element is burning
        end
        
        % Geometry update
        if FV_ignited(1,i) == 0
            Dint_p_t(2,i) = Dint_p_t(1,i);
            FV_dx(2,i) = FV_dx(1,i);
            FV_area(2,i) = FV_area(1,i);                  
        else
            Dint_p_t(2,i) = Dint_p_t(1,i) + 2*rdot_el*SIM(3);
            if FV_type(i) == 0 % End element
                FV_dx(2,i) = FV_dx(1,i) - (PROP(6)-1)*rdot_el*SIM(3);
                if FV_dx(2,i) < 0
                    FV_type(i) = -1; % Consumed element
                    % This code section is a workaround, but finds if the
                    % consumed end element is the first of final element 
                    % of a grain and changes the element type of the next
                    % or previous grain to an end grain and the current grain
                    % to consumed. The main difficult is to be able to work
                    % for multiple BATES grains. Even with this vast switch 
                    % makes it compatible only up to 4 BATES grain.
                    
                    switch i % End element is the first of a grain
                        case k1(1)
                            k1(1)=k1(1)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case SIM(4)/PROP(4)+k1(2)
                            k1(2)=k1(2)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 2*SIM(4)/PROP(4)+k1(3)
                            k1(3)=k1(3)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 3*SIM(4)/PROP(4)+k1(4)
                            k1(4)=k1(4)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 4*SIM(4)/PROP(4)+k1(5)
                             k1(5)=k1(5)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 5*SIM(4)/PROP(4)+k1(6)
                            k1(6)=k1(6)+1; FV_type(i+1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        otherwise
                    end
                    
                    switch i % End element is the last of a grain
                        case SIM(4)/PROP(4)-k2(1)
                            k2(1)=k2(1)+1; FV_type(i-1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 2*SIM(4)/PROP(4)-k2(2)
                             k2(2)=k2(2)+1; FV_type(i-1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 3*SIM(4)/PROP(4)-k2(3)
                            k2(3)=k2(3)+1; FV_type(i-1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 4*SIM(4)/PROP(4)-k2(4)
                            k2(4)=k2(4)+1; FV_type(i-1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        case 5*SIM(4)/PROP(4)-k2(5)
                            k2(5)=k2(5)+1; FV_type(i-1) = 0;
                            mdot_gen_el(i) = 0; Qw_el(i) = 0;
                        otherwise
                    end
                    
                else
                end
                FV_area(2,i) = pi*Dint_p_t(2,i)*FV_dx(2,i) - pi*Dint_p_t(2,i)*rdot_el*SIM(3) -...
                pi/4*(Dint_p_t(2,i)^2-Dint_p_t(2,i)^2); 
                % An end element losses area due to burning from the ends,
                % irrelevant for low N
            elseif FV_type(i) == -1  % Consumed element
                % Is left unused
            else
                FV_dx(2,i) = FV_dx(1,i);
                FV_area(2,i) = pi*Dint_p_t(2,i)*FV_dx(2,i);           
            end  
        end 
        
    end % XXXXX END OF ELEMENTWISE XXXXX
     
    Q_tot = sum(Qw_el); % Summation of total lost heat
    mdot_gen = sum(mdot_gen_el); % Summation of total generated mass flow
    
    if 1.001*(mdot_gen+mdot_ign) >= choked_nozzle(Dt_no_t(it),P_ch(it),T_ch(it),SIM(2),PROP(15),PROP(26),1)
        % Factor 1.001 is a workaround to avoid too close choked condition
        % Choked nozzle
        mdot_noz = choked_nozzle(Dt_no_t(it),P_ch(it),T_ch(it),SIM(2),PROP(15),PROP(26),1);
    else
        % Unchoked nozzle
        mdot_noz = choked_nozzle(Dt_no_t(it),P_ch(it),T_ch(it),SIM(2),PROP(15),PROP(26),0);
    end
    
    dT_dt = PROP(26)*T_ch(it)/(P_ch(it)*V_free*PROP(14))*(mdot_gen*PROP(15)*PROP(14)*T_flame+...
        mdot_ign*IGN(1,3)*IGN(1,4)*IGN(1,2)-mdot_noz*PROP(15)*PROP(14)*T_ch(it)-T_ch(it)*(PROP(14)*mdot_gen+...
        IGN(1,4)*mdot_ign-PROP(14)*mdot_noz)) + Q_tot/PROP(14);
    dP_dt = PROP(26)*T_ch(it)/V_free*(mdot_gen+mdot_ign-mdot_noz)+P_ch(it)/T_ch(it)*dT_dt-P_ch(it)*mdot_gen/(V_free*PROP(9));
  
    T_ch(it+1) = T_ch(it)+dT_dt*SIM(3);
    P_ch(it+1) = P_ch(it)+dP_dt*SIM(3);
    % Euler method for both chamber pressure and temperature   
    
    if T_ch(it+1) > T_flame
        T_ch(it+1) = T_flame; % Chamber cannot be hotter than the flame
    end 
    
    Dt_no_t(it+1) = Dt_no_t(it) + 0.7893*PROP(23)*(P_ch(it+1)/10^9)^(0.8)*10^(-3)*SIM(3) + dDt_no*SIM(3); 
    % Nozzle diameter erosion due to model + experimental 6% diameter
    % increase
  
    Pe = P_ch(it)/PR;   
    CF(it) = eff_tot*sqrt((2*PROP(15)^2)/(PROP(15)-1)*(2/(PROP(15)+1))^((PROP(15)+1)/(PROP(15)-1))*(1-(1/PR)^((PROP(15)-1)/PROP(15)))) + (Pe-SIM(2))*NOZ(3)^2/(P_ch(it)*Dt_no_t(it)^2);
    if CF(it) < 0
        CF(it) = 0;
    end
    if mdot_gen > 0
        cstar(it) = PROP(22)*P_ch(it)*pi/4*Dt_no_t(it)^2/(mdot_gen);
        % Chamber property
    else
        cstar(it) = 0;
    end
    if cstar(it) > 10^10
        cstar(it) = 0;
    end
    
    T(it) = CF(it)*P_ch(it)*pi/4*Dt_no_t(it)^2;
    
    Isp(it) = cstar(it)*CF(it)/9.81;
    
    % Provides output to the user each 10 ms for the first 200 ms and each
    % 100 ms after that
    if t_vec(it) <= 0.2
        if abs(t_vec(it)-0.01*(it1-1)) < 10^(-10)
            fprintf('Results at t = %4.2f s / ',t_vec(it));
            fprintf('Chamber pressure = %4.2f bar / ',P_ch(it)*10^(-5));
            fprintf('Chamber temperature = %4.2f K / ',T_ch(it));
            fprintf('Propellant mass left = %4.2f kg',round(m_p,3));
            fprintf('\n'); % New line
            it1 = it1 + 1;
        end
    elseif t_vec(it) > 0.2
        if abs(t_vec(it)-0.1*(it2-1)) < 10^(-10)
            fprintf('Results at t = %4.2f s / ',t_vec(it));
            fprintf('Chamber pressure = %4.2f bar / ',P_ch(it)*10^(-5));
            fprintf('Chamber temperature = %4.2f K / ',T_ch(it));
            fprintf('Propellant mass left = %4.2f kg',round(m_p,3));      
            fprintf('\n');
            it2 = it2 + 1;
        end      
    else
    end
    
    if P_ch(it) > stress_interpolation(SIM(1),yield_ch_T,yield_ch,uts_ch_T,uts_ch)*(CH(1)-CH(2))/CH(1)
        disp('----------------------------------');
        disp('Failure of the combustion chamber.')
        f_fail = true;
        break;
    elseif P_ch(it)*pi/4*CH(2)^2 > CH(4)
        disp('----------------------------------');
        disp('Failure of the bulkhead.')
        f_fail = true;
        break;
    elseif P_ch(it)*pi/4*(CH(2)^2-Dt_no_t(it)^2) > NOZ(6) || T_ch(it) > NOZ(10)
        disp('----------------------------------');
        disp('Failure of the nozzle.')
        f_fail = true;
        break;
    elseif T_ch(it) > NOZ(10)
        disp('----------------------------------');
        disp('Failure of the nozzle (temperature).')
        f_fail = true;
    else   
    end
    
    if (m_p < 10^(-6)) || (any(PROP(1)-Dint_p_t(2,:)<0))
    % Stops simulation for consumed mass or elements regression entering
    % the negative
        it_off = it;
        break;
    end
    
    it=it+1;
end

%% CALC: Tail-off
while (P_ch(it-1) > SIM(2)) && ~f_fail
    P_ch(it) = P_ch(it_off)*exp(-(PROP(26)*T_ch(it_off)*pi/4*Dt_no_t(it_off)^2)/(V_free*cstar(it_off))*SIM(3)*(it-it_off));
    % Nozzle erosion is negligible here
    Pe(it) = P_ch(it)/PR;
    CF(it) = eff_tot*sqrt((2*PROP(15)^2)/(PROP(15)-1)*(2/(PROP(15)+1))^((PROP(15)+1)/(PROP(15)-1))*(1-(1/PR)^((PROP(15)-1)/PROP(15)))) + (Pe(it)-SIM(2))*NOZ(3)^2/(P_ch(it)*Dt_no_t(it_off)^2);
    T(it) = eff_tot*CF(it)*P_ch(it)*pi/4*Dt_no_t(it_off)^2;
    Isp(it) = cstar(it_off)*CF(it)/9.81;
    if T(it) < 0 || CF(it) < 0
        CF(it) = 0;
        T(it) = 0;
    end
    
    it = it+1;
end
disp('-------------------------------------------------------------------------')
fprintf('Run is finished\n')
toc; % Execution time chronometer

%% OUTPUTS
if ~f_fail
cont = 0;
for j=1:it-1
    if it <= 50000 % Saves all results for low amount of data points
        results(2,cont+2) = t_vec(j);
        results(3,cont+2) = P_ch(j);
        results(4,cont+2) = T(j);
        results(5,cont+2) = CF(j);
        cont = cont + 1;        
    elseif t_vec(j) <= 0.2
        if rem(j,10) == 0 % Saved 1/10 data
            results(2,cont+2) = t_vec(j+1);
            results(3,cont+2) = P_ch(j+1);
            results(4,cont+2) = T(j+1);
            results(5,cont+2) = CF(j+1);
            cont = cont + 1;             
        elseif j == 1
            results(2,cont+2) = t_vec(j);
            results(3,cont+2) = P_ch(j);
            results(4,cont+2) = T(j);
            results(5,cont+2) = CF(j);
            cont = cont + 1;    
        else
        end               
    else
        if rem(j,100) == 0 % Saves 1/100 data
            results(2,cont+2) = t_vec(j+1);
            results(3,cont+2) = P_ch(j+1);
            results(4,cont+2) = T(j+1);
            results(5,cont+2) = CF(j+1);
            cont = cont + 1;             
        else
        end
    end
    
end
results(1,:) = zeros(1,length(results(2,:)));
results(1,2) = length(results(1,:))-1;
results(6,:) = zeros(1,length(results(2,:)));
results(6,2) = sum(T(1:it)*SIM(3));

Isp_av = results(6,2)/(9.81*pi/4*(PROP(1)^2-PROP(2)^2)*PROP(3)*PROP(4)*PROP(9));

disp('-------------------------------------------------------------------------')
disp('In which format do you want to save the results? (either .csv or .mat)');
file_type = input('.mat is a native Matlab format while .csv is more compatible: ','s');
if strcmpi('mat',file_type) || strcmpi('.mat',file_type) || strcmpi('m',file_type) || strcmpi('ma',file_type) % Just in case
    save('output.mat','results');
else % As default csv will be used
    filename = 'output.csv';
    writematrix(results,filename);
end
Post_processing(results); % Calls Post_processing with the results vector
else
end

%% Functions

% Intepolation of stress as function of temperature
function [yield_intpl,uts_intpl] = stress_interpolation(T,yield_T,yield,uts_T,uts) % Finds the uts and yield at a certain temperature
    if length(yield) == 1 % Only one tabulated value
        yield_intpl = yield;
        uts_intpl = uts;
    else
        yield_test = yield_T > T; 
        int_yield = find(yield_test,1);
        if yield_test == ones(1,length(yield_test)) % All data points are over the requested T
            yield_intpl = (yield(2)-yield(1))/(yield_T(2)-yield_T(1))*(T-yield_T(1))+yield(1);
            % The two first data points are used to create a line
        elseif yield_test == zeros(1,length(yield_test)) % All data points are below the requested T
            yield_intpl = (yield(end)-yield(end-1))/(yield_T(end)-yield_T(end-1))*(T-yield_T(end-1))+yield(end-1);
            % The two last data points are used to create a line
        else
            yield_intpl = (yield(int_yield)-yield(int_yield-1))/(yield_T(int_yield)-yield_T(int_yield-1))*(T-yield_T(int_yield-1))+yield(int_yield-1);
            % Normal linear interpolation
        end
        
        % The same is done for the ultimate tensile strength
        uts_test = uts_T > T; int_uts = find(uts_test,1);
        if uts_test == ones(1,length(uts_test))
            uts_intpl = (uts(2)-uts(1))/(uts_T(2)-uts_T(1))*(T-uts_T(1))+uts(1);
        elseif uts_test == zeros(1,length(uts_test))
            uts_intpl = (uts(end)-uts(end-1))/(uts_T(end)-uts_T(end-1))*(T-uts_T(end-1))+uts(end-1);
        else
            uts_intpl = (uts(int_uts)-uts(int_uts-1))/(uts_T(int_uts)-uts_T(int_uts-1))*(T-uts_T(int_uts-1))+uts(int_uts-1);
        end
    end
end

% Reduction of flame temperature due to combustion efficiency
function T_flame = flame_temp(T_flame_p,eff_comb)
    T_flame = T_flame_p*eff_comb;
    % A simple model to reduce adiabtic flame temperature is used
end

% Mass flows out of nozzle for choked and unchoked conditions
function mdot = choked_nozzle(Dt_no,P_ch,T_ch,Pamb,k_g,r_g,f_choked)
    At = pi/4*Dt_no^2;
    if f_choked == 1
        mdot = At*P_ch*sqrt(k_g/(r_g*T_ch)*(2/(k_g+1))^((k_g+1)/(k_g-1))); % we consider choked nozzle  
    else
        mdot = 0.85*At*sqrt(2*P_ch/(r_g*T_ch)*(P_ch-Pamb)); % we consider quadrant edge orifice
    end
    if ~isreal(mdot) || mdot < 0
        mdot = 0;
    end
end

% Erosive burning rate calculation
function rdot_ero = erosive_burning(dens_g,v_g,Dh,a_rdot,b_rdot,rdot_nom,dens_p) % Everything is with SI units
    a_rdot = a_rdot*(rdot_nom*Dh)^0.2/dens_p; % Correction given by Hein1
    G = dens_g*v_g;
    z = (b_rdot*dens_p*rdot_nom)/G;   
    rdot_ero = a_rdot*G^(0.8)/(Dh^0.2*exp(z)); % Lenoir-Robillard
end

% 4th Order Runge-Kutta for Twall differential equation
function twall_next = Twall_4th_order_RK (diff_p,h_c,h_r,lambda_p,Tamb,Twall0,Tg1,Tg0,dt_sim)
    % x --> t // y --> Twall
    l = lambda_p; d = diff_p; h_t = h_r + h_c; % Short forms for clarity
    h = dt_sim/100; % Step-size a hundreth. A tenth already yielded precise results
    x = 0:h:dt_sim; y = zeros(1,length(x));
    % Tg linearization
    m = (Tg1-Tg0)/dt_sim; n = Tg0;
    
    x(1) = 0; y(1) = Twall0;
    y_dot =@(x,y)(4*d*(h_t/l*(m*x+n-y))^3/(3*(2*h_t/l*(m*x+n-y)*(y-Tamb)+h_t/l*(m*x+n-Tamb)^2))); % HEIN REDUCED EXPRESSION
    for i=1:length(x)-1
        k1 = y_dot(x(i),y(i));
        k2 = y_dot(x(i)+.5*h,y(i)+.5*k1*h);
        k3 = y_dot(x(i)+.5*h,y(i)+.5*k2*h);
        k4 = y_dot(x(i)+h,y(i)+k3*h);
        y(i+1) = y(i)+((k1+2*k2+2*k3+k4)/6)*h;
    end
    twall_next = y(end);
end

% Finite volume discretization of the grain
function [FV_type,FV_length,FV_area] = FV_discr (L_ch,L0_p,N_p,N_sim,Dext_p,Dint_p)
    N0_p = ceil(N_sim/N_p); % Rounds up and assigns the number of FV that correspond to each grain
    %N_sim_old = N_sim; % Saves the initial amount of FVs
    N_sim = N0_p*N_p; % Corrects the number of FV finally used (if necesary, N_sim might be equal to N_sim_old)
    dL = L0_p/N0_p;
    dL_free = (L_ch-N_p*L0_p)/(2*N_p);
    FV_type_N = ones(1,N0_p); % N0_p
    FV_type_N(1) = 0; FV_type_N(end) = 0; % End elements
    FV_type = FV_type_N;
    
    if N_p > 1
        for i=1:N_p-1
            FV_type = [FV_type FV_type_N]; % Creates final vector of element type
        end
    else    
    end

    FV_length_N = zeros(N0_p,2); % Position at which the FV starts (both core and end elements) 
    for i=1:N0_p
        FV_length_N(i,:) = [dL_free+(i-1)*dL, dL_free+i*dL];
    end
    
    FV_length = FV_length_N;
    if N_p > 1 % Concatanted the FV_length
        for i=1:N_p-1
            FV_length = [FV_length; dL_free+FV_length(end,2)+FV_length_N];
        end
    else    
    end   
       
    FV_area = zeros(1,length(N_sim));
    for i = 1:N_sim
        FV_area(i) = pi*Dint_p*dL;
    end
end