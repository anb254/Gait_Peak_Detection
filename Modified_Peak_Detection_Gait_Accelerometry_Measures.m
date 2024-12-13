%This code uses stride segmented data to calculate gait features of
%interest. The input is a .mat file for which heel strikes have been
%identified. The output is an excel file that contains step time, HR, and
%magnitude of acceleration within each stride.

%Edited by Anna Bailes 10/31/23

%% LOAD FILES

close all; clear all; clc

for j=1
    subject= %Insert subject number
    fileData = strcat('P:\Gait and LBP\Aim 1\Data Processing\Raw Participant Data\000',num2str(subject),'\segmented_000',num2str(subject),'_2MWT Data.mat');
    data=load(fileData);

%% ASSIGN DATA VARIABLES AND FILTER DATA

    ap_acc=data.ap_acc; %z column, +ve=forward (changed sign in peak detection code)
    ml_acc=data.ml_acc;  %y column, +ve=right
    v_acc=data.v_acc; %x column, +ve=upward (changed sign in peak detection code) 
    loc=data.loc;
    time=data.time;
    
    %Apply low pass Butteworth filter 
    fs=62.5;
    fc=15; %almost all signal power is captured under 15 Hz (Kavanagh 2005)
    [b,a]=butter(4,fc/(fs/2)); %4th order LP filter, with cutoff freq normalized to nyquist freq (fs/2). Using cutoff freq of 15 (Kavanagh 2005)
    ap_acc_filt = filtfilt(b,a,ap_acc);
    ml_acc_filt = filtfilt(b,a,ml_acc);
    v_acc_filt = filtfilt(b,a,v_acc);
    

%% DEROTATION STEP (TO ALIGN WITH GRAVITY)

    acc_mat=horzcat(v_acc_filt,ml_acc_filt,ap_acc_filt); 
    acc_mean=mean(acc_mat);
    
    %XZ(rotation in sagittal plane)
    theta=atan(acc_mean(3)/acc_mean(1)); %using trig and geometry
    thetadeg=theta*(180/pi);
    Ry=[cos(theta),sin(theta);-sin(theta),cos(theta)]; %rotation matrix around y axis

    %XY (rotation in frontal plane)

    phi=atan(-acc_mean(2)/acc_mean(1)); %using trig and geometry
    phideg=phi*(180/pi);
    Rz=[cos(phi),-sin(phi);sin(phi),cos(phi)]; %rotation matrix around z axis

    R=[Ry(1,1),0,Ry(1,2);0,1,0;Ry(2,1),0,Ry(2,2)]*[Rz(1,1),Rz(1,2),0;Rz(2,1),Rz(2,2),0;0,0,1];
    acc_rot=R*acc_mat';
    acc_rot=acc_rot'; %this is the derotated data
    ap_acc_rot=acc_rot(:,3);
    ml_acc_rot=acc_rot(:,2);
    v_acc_rot=acc_rot(:,1);
    
%% MODIFY PEAK DETECTION METHOD

%Peak_Detection code only captured pendular movement of gait.
%No we used each of these indeces (loc) to find actual peak of acceleration
%signal, which is indicative of heel strike

%Find minima of AP acc signal
min_vect=[];
min_idx_vect=[];
c=1;
for i=1:length(loc)-1
    window=-ap_acc_rot(loc(i):loc(i+1)); %invert to get minima
    [m,idx]=max(window);
    idx=idx+loc(i)-1;
    min_vect(c,1)=m;
    min_idx_vect(c,1)=idx;
    c=c+1;
end

peak_vect=[];
peak_idx_vect=[];
c=1;
idx=[];


%Find all peaks between peak of 2 Hz filter waveform (loc) and acc signal
%minima (min_idx_vect)
for i=1:length(loc)-1
    if min_idx_vect(i)-loc(i)<3
        window=ap_acc_rot(loc(i-1):min_idx_vect(i+1));
    else
        window=ap_acc_rot(loc(i):min_idx_vect(i));
    end
    [p,idx]=findpeaks(window,'MinPeakHeight',0.02); %find peaks above 0.02
    if size(p,1)>1 %if multiple peaks, choose last one
        p=p(end);
        idx=idx(end);
    elseif isempty(p)
        [p,idx]=findpeaks(window,'MinPeakHeight',0.01); %if no big peaks, find peaks above 0.01
            if size(p,1)>1 %if multiple peaks, choose last one
                p=p(end);
                idx=idx(end);
            elseif isempty(p) % if still no peaks, the peak of windowed signal is the start (ie peak from 2 Hz filter)
                p=ap_acc_rot(loc(i));
                idx=1;
            end
    end
    idx=idx+loc(i)-1;
    peak_vect(c,1)=p;
    peak_idx_vect(c,1)=idx;
    c=c+1;
end 

folder=strcat('P:\Gait and LBP\Aim 1\Data Processing\Processed Data\Plots\000',num2str(subject));
mkdir(folder);

figure();
plot(ap_acc_rot,'b')
hold on
plot(loc,ap_acc_rot(loc),'k*')
plot(peak_idx_vect,ap_acc_rot(peak_idx_vect),'ro')
plot(min_idx_vect,ap_acc_rot(min_idx_vect),'b*')
title(subject)
fig_name=num2str(subject);
saveas(gcf,strcat(folder,'/',fig_name),'png');
saveas(gcf,strcat(folder,'/',fig_name),'fig');
% close

%% CALCULATE GAIT METRICS

% STEP TIME

    timepks = time(peak_idx_vect); %time vector is in seconds
    steptime = diff(timepks); %each peak corresponds to one step

    %Calculate step time for each side (left and right)
    if rem(length(steptime),2)==1 %this means there were an ODD number of steps identified
        steptime0 = steptime(1:2:end-2); %removing last stride so that there are an equal number of steps per side
        steptime1=steptime(2:2:end-1);
    elseif rem(length(steptime),2)==0 %this means there were an EVEN number of steps identified
        steptime0 = steptime(1:2:end-1);
        steptime1=steptime(2:2:end);
    end

% HARMONIC RATIO
 
    stride_loc=peak_idx_vect(1:2:end); %every other peak is used to identify a single stride

    %Call StrideHarmRatio function to calculate harmonic ratios for each stride

    HR_AP_vect=[];
    HR_V_vect=[];
    HR_ML_vect=[];

    %HR AP
    for i=1:length(stride_loc)-1
        stride=ap_acc(stride_loc(i):stride_loc(i+1));
        HR_vect_new=StrideHarmRatio(stride,1);
        HR_AP_vect=[HR_AP_vect,HR_vect_new];
    end
    HR_AP_vect=HR_AP_vect';

    %HR ML
    for i=1:length(stride_loc)-1
        stride=ml_acc(stride_loc(i):stride_loc(i+1));
        HR_vect_new=StrideHarmRatio(stride,0);
        HR_ML_vect=[HR_ML_vect,HR_vect_new];
    end

    HR_ML_vect=HR_ML_vect';

    %HR V
    for i=1:length(stride_loc)-1
        stride=v_acc(stride_loc(i):stride_loc(i+1));
        HR_vect_new=StrideHarmRatio(stride,1);
        HR_V_vect=[HR_V_vect,HR_vect_new];
    end

    HR_V_vect= HR_V_vect';

%% SAVE VARIABLES

    len=length(steptime0);
    subj_v=subject*ones(len,1);

    fileName=strcat('P:\Gait and LBP\Aim 1\Data Processing\Processed Data\00',num2str(subject),'_Gait_metrics.xls');
  
    header={'Subject','Step Time 0', 'Step Time 1','HR_AP','HR_ML','HR_V'};
    mat=[subj_v,steptime0,steptime1,HR_AP_vect,HR_ML_vect,HR_V_vect];
    writecell(header,fileName);
    writematrix(mat,fileName,'WriteMode','append');
end
