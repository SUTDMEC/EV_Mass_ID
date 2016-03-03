function out = RT_PI_ImIEV(in)
% try

% add 'in.debug' to in structure to plot and debug various options 
if ~any(strcmp('debug',fieldnames(in)))
    in.debug=0;
end

if ~any(strcmp('mass_guess',fieldnames(in)))
    true_mass= (1120 + 85) * 1.2885; %1.2885 is factor for apparent mass (only Erik) kerb weight = 1110 (measured @ 1120 on 17 apr 2013)
    % true_mass=(1256 + 140 + 15 + 10) * 1.2885; %1.2885 is factor for apparent mass (whole fam + silvie + gear)
else
    true_mass=(1120+ in.mass_guess) * 1.2885;
end


eff=0.8; %naive assumed motor efficiency

tau = smooth(in.tau,30); % Moving average smoothing
dtaudtime = diff(tau)./diff(in.tautime); %first derivative of torque (in this case power)
dtaudtime=[99 dtaudtime'];
dtaudtime(tau==0)=99; %mark elements where torque is 0 for removal
x_ind=zeros(size(dtaudtime));

%NOTE: TAU = POWER (kW) NOT TORQUE (N-m)

%power derivative limit sets how stable to noise the signal is
tau_d_limit=3.5;

x_ind((tau_d_limit>dtaudtime&dtaudtime>-tau_d_limit))=10; %power first derivative < 1 - i.e. signal is stable, slope = 0, set = flag of 10

%find consecutive values 
thresh=10; % threshold for # of consecutive stable counts before stable zone is found
startIndex = strfind(x_ind, 10.*ones(1,thresh)); % count the flags of 10
if ~isempty(startIndex) % see if a stable zone is found
    indices = unique( bsxfun(@plus, startIndex', 0:thresh-1) )';
    iplot=11.*ones(size(indices));
else
    out=in;
    out.success=0;
    return
end

%compare with velocity to find acceleration events
vel=smooth(in.vel./3.6,5); % Moving average smoothing
acc=diff(vel)./diff(in.veltime); %acceleration
acc= [0 acc'];
v_ind=zeros(size(acc));

%power derivative limit sets how stable to noise the signal is
v_d_limit_high=3;
v_d_limit_low=0.5;
v_ind((v_d_limit_high>acc&acc>v_d_limit_low))=20; %acceleration flags of 20 between the predefined limits

% keyboard
%find consecutive values for acceleration
vthresh=15; % threshold for # of consecutive stable counts before stable zone is found
vstartIndex = strfind(v_ind, 20.*ones(1,vthresh)); %count the flags of 20
if (~isempty(vstartIndex))
    vindices = unique( bsxfun(@plus, vstartIndex', 0:vthresh-1) )';
    viplot=15.*ones(size(vindices));
else
    out=in;
    out.success=0;
    return
end

%find indicies which overlap

rounded_tau=round(in.tautime(indices).*10)./10;
rounded_vel=round(in.veltime(vindices).*10)./10;

comm_time=intersect(rounded_tau,rounded_vel);
commpl=18.*ones(size(comm_time));
% keyboard

%% This method uses average torque/accel - depreciated
%find the places where there are sets of complimentary times
% clust=find(diff(comm_time)>10);
% 
% %run through the clusters to calculate mass at each one
% comm_start=1;
% for i=1:length(clust)+1
%     % Find average power/torque and acceleration for each cluster of 'good'
%     % points
%     if i<length(clust)
%         comm_end=clust(i);
%     else
%         comm_end=length(comm_time);
%     end
%     
%     [~, indstart_tau]=min(abs(in.tautime-comm_time(comm_start)));
%     [~, indend_tau]=min(abs(in.tautime-comm_time(comm_end)));
% 
%     [~, indstart_vel]=min(abs(in.veltime-comm_time(comm_start)));
%     [~, indend_vel]=min(abs(in.veltime-comm_time(comm_end)));
%     
%     indmid=floor((comm_start+comm_end)/2);
%     ave_tau(i)=mean(tau(indstart_tau:indend_tau)); %actually average power in this case, in kW
%     ave_acc(i)=mean(acc(indstart_vel:indend_vel)); %average acceleration in m/s^2
%     mid_time(i)=comm_time(indmid);
%     if i<length(clust)
%         comm_start=clust(i+1);
%     else
%         comm_start=clust(i-1);
%     end
%     
%     % calculate mass based on newton's second law
%     v_ang(i)= (INCOMPLETE)
%     torque(i)=ave_tau(i)*1000/v_ang(i);
%     
% end

%% Use instantaneous speed/power/acceleration to do the calculation at the 'common steady state points'

for j=1:length(comm_time)
    [~, ind_vel(j)]=min(abs(in.veltime-comm_time(j)));
    [~, ind_pow(j)]=min(abs(in.tautime-comm_time(j)));
    %ind_vel(j)=find(comm_time(j)==in.veltime);
    %ind_pow(j)=find(comm_time(j)==in.tautime);
end

% [~, ind_vel]=min(abs(in.veltime-comm_time));
% [~, ind_pow]=min(abs(in.tautime-comm_time));
if ~exist('ind_vel','var')
    out.text='no data';
    out.success=0;
    return
end
v_inst=vel(ind_vel);
p_inst=tau(ind_pow).*eff;
a_inst=acc(ind_vel);

gear_red=7.065;
tire_diameter=0.5695;

rpm_inst=v_inst./tire_diameter.*gear_red./2./pi().*60;
rps_inst=rpm_inst./60;
ang_vel=rps_inst.*2.*pi();
torq=p_inst.*1000./ang_vel;
force=torq.*gear_red./tire_diameter;
mass=force./a_inst';

% handle physical realizability 
ind_acc_vio=find(a_inst>2.7 | a_inst<0); %acceleration less that 10 s 0-100 kph (0 - 27.7 m/s)

[M,N] = size(ind_acc_vio);
if M<N
    ind_acc_vio=ind_acc_vio';
end

ind_p_vio=find(p_inst>60 | p_inst<0); %power > 47 kW, less than 0
ind_tau_vio=find(torq>220 | torq<0); %torque > 180 N-m, less than 0
ind_mass_vio=find(mass<1450 | mass > 1800); %1400 , 1800 is 50 kg > the rated GVM multiplied by apparent mass factor http://www.mitsubishi-cars.co.uk/imiev/specifications.aspx

out.violations=[length(ind_acc_vio);length(ind_p_vio);length(ind_tau_vio);length(ind_mass_vio)];

out.mass_not_parse=mass;
combined=[ind_acc_vio;ind_p_vio;ind_tau_vio;ind_mass_vio];
mass_parse=mass;
mass_parse(combined)=[];
parse_time=comm_time;
parse_time(combined)=[];


out.vel=vel;
out.mass=mass_parse;
out.parse_time=parse_time;
out.force=force;
out.torq=torq;
out.comm_time=comm_time;
error = (mass_parse - true_mass)./true_mass;

out.error = error;
out.true_mass=true_mass;

%handle empty mass guess case
if ~isempty(mass_parse)
    out.cum_mass_parse=cummean(mass_parse);
    out.cum_error=cummean(error);
else
    out.cum_mass_parse=0;
    out.cum_error=0;
end

%%APPLY THIS METHOD TO ESTIMATE ERROR
%http://www.mathworks.com/help/dsp/ref/msepred.html
%%

if ~isempty(out.mass) && ~isempty(out.error)
    out.success=1;
else
    out.success=0;
end
%% Plot results
% keyboard


%Truncated results plot (mass and error only)
if in.debug
    figure
    %full results plotting
    if ~isempty(iplot)
        plot(in.tautime,tau,'k','LineWidth',1.5);
        hold on
        plot(in.veltime,vel,'LineWidth',1.5);
        %plot(in.tautime,dtaudtime,'k'); 
        %plot(in.tautime,x_ind,'r'); 
        %plot(in.tautime(indices),iplot,'*k'); 
        %plot(in.veltime(vindices),viplot,'*g');
        %plot(comm_time,commpl,'*m');
    %     plot(mid_time,ave_tau,'sr','Markersize',12);
    %     plot(mid_time,ave_acc,'sb','Markersize',12);

    % figure
    plot(parse_time,mass_parse./100,'sr','Markersize',12);
    plot(comm_time,force./100,'sb','Markersize',12);
    plot(comm_time,a_inst,'sk','Markersize',12);
    plot(comm_time,torq./10,'sg','Markersize',12);
    % plot(comm_time,mass,'sb','Markersize',12, 'MarkerFaceColor','b');
    plot(parse_time,error.*100,'sr','Markersize',12, 'MarkerFaceColor','r');
    legend({'tau (N-m)' 'spd (km/h)' 'mass/100 (kg)' 'force/100 (N)' 'accel (m/s^2)' 'torq/10 (N-m)' 'error (%)' });
    title('Time-series Data');
    ylabel('Scaled values');
    xlabel('Time (s)');
    set(gca,'LineWidth',1.5)
    figure; 
    plot(out.mass,'*');
    hold
    plot(out.cum_mass_parse,'r','LineWidth',1.5);
    plot(1:length(out.mass),out.true_mass,'*k');
    title('Identified Mass');
    ylabel('Mass (kg)');
    xlabel('Event');
    legend({'Ident. Mass' 'Cum. ave.' 'True mass' });
    set(gca,'LineWidth',1.5)
    ylim([0 2000])
    end
elseif in.plot_light
    plot(parse_time,error.*100,'sr','Markersize',12, 'MarkerFaceColor','r');
end
% keyboard
% catch
%     a=lasterror;
%     keyboard;
% end