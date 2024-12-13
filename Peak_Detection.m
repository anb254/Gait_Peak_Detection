%This code requires an input of 3D linear trunk accelerations and uses a
%"peak detection" stride segmentation method to separate continuous gait
%data into individual steps. The output is a .mat file that
%contains time, 3D linear acceleration, and time indexes corresponding to
%stide segmentation (i.e. heel strike) events

%Method modified from Zijlstra 2003 (peak detection method)

%X=superior->inferior (vertical)
%Y=left->right (ML)
%Z=anterior->posterior (AP)

%Last edited by Anna Bailes 1/3/2023


%% IDENTIFY HEEL STRIKES

%Use AP acc to identify peaks that correspond to heel strike 

%Apply low pass filter to capture pendular movement
[b,a]=butter(4,2/(fs/2)); %4th order LP filter, with cutoff freq normalized to nyquist freq (fs/2). Using cutoff freq of 2 to capture pendular movement of gait
ap_acc_filt = filtfilt(b,a,ap_acc); %use anterior-posterior trunk accelerations

[pks, loc] = findpeaks(ap_acc_filt);

%% REMOVE ACCELERATION/DECELERATION PHASES OF 2MWT

%Remove first/last 5 steps of trial for acceleration/deceleration phase.

figure(3)
hold on
plot(ap_acc,'r') %unfiltered AP acc
plot(ap_acc_filt,'b')  %filtered AP acc
plot(loc,ap_acc_filt(loc),'k*')
set(gcf, 'WindowState', 'maximized');
legend('Unfiltered','Filtered','Steps')

%Zoom in to identify first "true step" 
zoom on;
pause()
zoom off;

[endpt1,~]=ginput(1);
zoom out;

%Zoom in to identify last "true step"
zoom on;
pause()
zoom off;

[endpt2,~]=ginput(1);
close all;
 
%Find closest asterisk to points identified from plot
[~,closestFirst]=min(abs(loc-endpt1));
[~,closestStop]=min(abs(loc-endpt2));

%Create new loc vector with acceleration and deceleration phases (5 steps) removed
loc=loc(closestFirst+5:closestStop-5);
step_count=length(loc); %number of steps


%% CHECK PLOT FOR MISSED STEPS

%Check new loc vector plot
h=figure(4);
hold on
plot(ap_acc,'r') %unfiltered AP acc
plot(ap_acc_filt,'b') %filtered AP acc
plot(loc,ap_acc_filt(loc),'k*')
set(gcf, 'WindowState', 'maximized');
legend('Unfiltered','Filtered','Steps')
title(strcat(subject,' Stride Segmented Accel.- Step Count = ',num2str(step_count)));
zoom on
pause()  

% Add any steps that were missed
mis=input('Were any steps missed? (Y/N) ','s');

if strcmp(mis,'N')
    disp('No steps missed')

elseif strcmp(mis,'Y')
    disp('Identify missed steps')
    h=figure;
    hold on
    plot(ap_acc,'r') %unfiltered AP acc
    plot(ap_acc_filt,'b') %filtered AP acc
    plot(loc,ap_acc_filt(loc),'k*')
    set(gcf, 'WindowState', 'maximized','currentchar',' ');
    legend('Unfiltered','Filtered','Steps')
    zoom on 
    pause()

    while ishandle(h)
        [new_step,~]=ginput(1); %Manually input missed steps
        new_step=round(new_step);
        loc(end+1)=new_step;
        pause()
    end 
 
     loc=sort(loc);
end

step_count_new=length(loc);

%% FINAL PLOT AND SAVE

figure(6)
hold on
plot(ap_acc,'r') %unfiltered AP acc
plot(ap_acc_filt,'b') %filtered AP acc
plot(loc,ap_acc_filt(loc),'k*')
set(gcf, 'WindowState', 'maximized');
legend('Unfiltered','Filtered','Steps')
title(strcat(subject,' Stride Segmented Accel.- Step Count = ',num2str(step_count_new)));

%Save plot to png and fig
fig_name='Stride Segmented Accel Plot';
saveas(gcf,strcat(folder,'/',fig_name),'png');
saveas(gcf,strcat(folder,'/',fig_name),'fig');

%Save acceleration vectors and location (index) vector to .mat file
savename=strcat(folder,'/','segmented_',fileData);
save(savename,'time','ap_acc','ml_acc','v_acc','loc');


