clear;
close all;
d = 1:5000;
a = 3;
b = 4;
for i = 1:5000
    SA(i) = 4*asin(1*b/sqrt((a*a+4*d(i)*d(i))*(b*b+4*d(i)*d(i))));
end
loglog(SA)
ylabel('Solid Angle')
xlabel('Distance')
title('solid angle of wall; size of wall=3*4')

%% Preprocess Transmittance Data
clear all;
close all;
load('transmittance_7_14')
% P(x=k)
temp_diff=10;
channels=100;
TIME=0.001;
area=5e-6^2;
distance=100;
expt = Expected_Photon(temp_diff,channels,TIME,area);

% Fisher information of 10 channels
wavelength = linspace(7,14,channels);
figure(1)
subplot(2,2,1)
result_10 = FisherInformation(temp_diff, 10,TIME,...
    area,transmittance);
plot(result_10)
title('Fisher information of 10 channels')
% Fisher information of 100 channels
subplot(2,2,2)
result_100 = FisherInformation(temp_diff, 100,TIME,...
    area,transmittance);
plot(result_100)
title('Fisher information of 100 channels')
% Fisher information of 500 channels
subplot(2,2,3)
result_1000 = FisherInformation(temp_diff, 1000,TIME,...
    area,transmittance);
plot(result_1000)
title('Fisher information of 1000 channels')

subplot(2,2,4)
plot(transmittance(:,1),1-transmittance(:,2))
title('1-transmittance')

% study the relationship between channels and total F-information
figure(2)
num_channel=1:100;
for i=1:100
    result(i) = sum(FisherInformation(temp_diff,num_channel(i),TIME,...
    area,transmittance));
end
plot(num_channel,result)
title('Relationship between number of channels and F-information')

% figure(3)
% result_4 = FisherInformation(temp_diff, 4,time,...
%     area,transmittance);
% subplot(1,2,1)
% plot(result_4)
% title('Fisher information of 4 channels')
% 
% result_5 = FisherInformation(temp_diff, 5,time,...
%     area,transmittance);
% subplot(1,2,2)
% plot(result_5)
% title('Fisher information of 5 channels')

% study the relationship between tem_diff and total F-informataion
figure(4)
diff=1:100;
num_channel = 100;
result = zeros(1,100);
for i=1:100
    result(i) = sum(FisherInformation(diff(i),num_channel,TIME,...
    area,transmittance));
end
plot(diff,result)
title('Relationship between temperature difference and F-information')


%% Find the best way to seperate to 2 bands
% close all;
load('transmittance_7_14')
temp_diff=10;
channels=2;
TIME=0.001;
area=5e-6^2;
distance=100;
T_obj = 300; 
T_atm = T_obj-10;
SUBCHANNELS = 1000;
wavelength = linspace(7,14,1000+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,SUBCHANNELS,TIME,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,SUBCHANNELS,TIME,area);
d = 100;
index = round(linspace(1,SUBCHANNELS+1,channels+1));


tem_result_up = zeros(1,1000);
tem_result_down = zeros(1,1000);

for i=1:1000
    c1 = ExpectedPhoton_atm(i);
    c2 = ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i);
    [~,trans_index] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index,2);
%         before change to sum Fisher information
%         tem_result(i) = c2^2*transmittance^(2*d)*(log(transmittance))^2/...
%         (c1+c2*transmittance^d); 


    tem_result_up(i) = c2*transm_single^d*log(transm_single);
    tem_result_down(i) = c1+c2*transm_single^d;
end
U = sum(tem_result_up);
D = sum(tem_result_down);
for ch = 1:channels
    tem_up = sum(tem_result_up(1,index(ch):index(ch+1)-1))^2;
    tem_down = sum(tem_result_down(1,index(ch):index(ch+1)-1));
    control_group(ch) = tem_up/tem_down;
end
result = zeros(2,SUBCHANNELS-1);
a = zeros(1,SUBCHANNELS-1);
b = zeros(1,SUBCHANNELS-1);
for i = 1:SUBCHANNELS-1
    tem_up = sum(tem_result_up(1,1:i))^2;
    tem_down = sum(tem_result_down(1,1:i));
    result(1,i) = tem_up/tem_down;
    a(1,i) = sum(tem_result_up(1,1:i))/U;
    b(1,i) = tem_down/D;
    
    tem_up = sum(tem_result_up(1,i+1:SUBCHANNELS))^2;
    tem_down = sum(tem_result_down(1,i+1:SUBCHANNELS));
    result(2,i) = tem_up/tem_down;
%     ratio_a(i) = D*b(i)/U;
%     ratio_a2(i) = D*b(i)/(2*b(i)-U);
%     ratio_b(i) = a(i)*U/D;
end
% tem_result_up = sum(tem_result_up)^2;
% tem_result_down = sum(tem_result_down);
% result = tem_result_up/tem_result_down;

figure()
hold on
plot(sum(result))
% plot([1,subchannels-1],[0.0046,0.0046],'r')
[~,biggest_index]=max(sum(result));
% plot([biggest_index,biggest_index],[0.004,0.01])
title('relationship between n and the overall F-information')
xlabel('Position of the interval')
ylabel('Overall F-information')
hold off
figure()
subplot(2,1,1)
plot(result(1,:))
title('the Fisher-information of the first channel')
subplot(2,1,2)
plot(result(2,:))
title('the Fisher-information of the second channel')
%% Find the best way to seperate to 3 bands
% close all;
load('transmittance_7_14')
temp_diff=10;
channels=3;
TIME=0.001;
area=5e-6^2;
T_obj = 300; 
T_atm = T_obj-10;
SUBCHANNELS = 1000;
wavelength = linspace(7,14,1000+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,SUBCHANNELS,TIME,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,SUBCHANNELS,TIME,area);
d = 100;
index = round(linspace(1,SUBCHANNELS+1,channels+1));


tem_result_up = zeros(1,1000);
tem_result_down = zeros(1,1000);

for i=1:1000
    c1 = ExpectedPhoton_atm(i);
    c1_record(i) = c1;
    c2 = ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i);
    c2_record(i) = c2;
    [~,trans_index] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index,2);
%         before change to sum Fisher information
%         tem_result(i) = c2^2*transmittance^(2*d)*(log(transmittance))^2/...
%         (c1+c2*transmittance^d); 
    transm_record(i) = transm_single;


    tem_result_up(i) = c2*transm_single^d*log(transm_single);
    tem_result_down(i) = c1+c2*transm_single^d;
end
U = sum(tem_result_up);
D = sum(tem_result_down);

% calculate F-infor of control group
% for ch = 1:channels
%     tem_up = sum(tem_result_up(1,index(ch):index(ch+1)-1))^2;
%     tem_down = sum(tem_result_down(1,index(ch):index(ch+1)-1));
%     control_group(ch) = tem_up/tem_down;
% end

result = zeros(3,498501);
index_ij = zeros(2,498501);
k=1;
for i = 2:999
    for j= i+1:1000
        tem_up = sum(tem_result_up(1,1:i-1))^2;
        tem_down = sum(tem_result_down(1,1:i-1));
        result(1,k) = tem_up/tem_down;

        tem_up = sum(tem_result_up(1,i:j-1))^2;
        tem_down = sum(tem_result_down(1,i:j-1));
        result(2,k) = tem_up/tem_down;
        
        tem_up = sum(tem_result_up(1,j:1000))^2;
        tem_down = sum(tem_result_down(1,j:1000));
        result(3,k) = tem_up/tem_down;
        
        index_ij(:,k) = [i;j];
        
        k = k+1;
    end
end
% tem_result_up = sum(tem_result_up)^2;
% tem_result_down = sum(tem_result_down);
% result = tem_result_up/tem_result_down;

figure()
hold on
plot(sum(result))
[~,biggest_index]=max(sum(result));
% plot([biggest_index,biggest_index],[0.004,0.01])
title('relationship between n and the overall F-information')
xlabel('Position of the interval')
ylabel('Overall F-information')
hold off


%% Find the best way to seperate to multiple bands
close all;
clear all;
load('transmittance_7_14')
temp_diff=10;
channels=50;
TIME=0.001;
area=5e-6^2;
distance=100;
T_obj = 300; 
T_atm = T_obj-10;
SUBCHANNELS = 1000;
wavelength = linspace(7,14,SUBCHANNELS+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,SUBCHANNELS,TIME,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,SUBCHANNELS,TIME,area);
d = 100;
index = round(linspace(1,SUBCHANNELS+1,channels+1));


tem_result_up = zeros(1,SUBCHANNELS);
tem_result_down = zeros(1,SUBCHANNELS);

for i=1:SUBCHANNELS
    c1 = ExpectedPhoton_atm(i);
    c2 = ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i);
    [~,trans_index] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index,2);
%         before change to sum Fisher information
%         tem_result(i) = c2^2*transmittance^(2*d)*(log(transmittance))^2/...
%         (c1+c2*transmittance^d); 


    tem_result_up(i) = c2*transm_single^d*log(transm_single);
    tem_result_down(i) = c1+c2*transm_single^d;
end
U = sum(tem_result_up);
D = sum(tem_result_down);

for ch = 1:channels
    tem_up = sum(tem_result_up(1,index(ch):index(ch+1)-1))^2;
    tem_down = sum(tem_result_down(1,index(ch):index(ch+1)-1));
    control_group(ch) = tem_up/tem_down;
end

result = zeros(2,SUBCHANNELS-1);
a = zeros(1,SUBCHANNELS-1);
b = zeros(1,SUBCHANNELS-1);

for ch = 2:30
    index = round(linspace(1,SUBCHANNELS+1,ch+1));
    for j = 2:(ch)
        ind = 1;
        subresult = zeros(2,1);
        for i = index(j-1)+1:index(j+1)-2
            tem_up = sum(tem_result_up(1,index(j-1)+1:i))^2;
            tem_down = sum(tem_result_down(1,index(j-1)+1:i));
            subresult(1,ind) = tem_up/tem_down;
    %         a(1,i) = sum(tem_result_up(1,1:i))/U;
    %         b(1,i) = tem_down/D;

            tem_up = sum(tem_result_up(1,i+1:index(j+1)-1))^2;
            tem_down = sum(tem_result_down(1,i+1:index(j+1)-1));
            subresult(2,ind) = tem_up/tem_down;
            ind = ind+1;
        end
        [~,biggest_index]=max(sum(subresult));
        index(j)=biggest_index+index(j-1);
    end
    optimized_result(ch) = sum(FisherInformation(temp_diff, ch,TIME,area,transmittance,index));
    index_1 = round(linspace(1,SUBCHANNELS,ch+1));
    referece(ch) = sum(FisherInformation(temp_diff, ch,TIME,area,transmittance,index_1));
end

figure()
plot(1:30,(optimized_result(find(optimized_result))-referece(find(optimized_result)))./referece(find(optimized_result)))
title('Improvement of Optimized Channels')

figure()
optimized_result(1) = 0.0045;
plot(1:30,(optimized_result))
load("Optimized")
hold on 
plot(1:30,(optimized_result))
legend('old method','new method')
title('Optimized Fisher-information of 1~30 Channels')
% 
% 
% figure()
% subplot(2,1,1)
% hold on
% plot(sum(result))
% plot([1,subchannels-1],[0.0046,0.0046],'r')
% [~,biggest_index]=max(sum(result));
% % plot([biggest_index,biggest_index],[0.004,0.01])
% title('relationship between n and the overall F-information')
% xlabel('Position of the interval')
% ylabel('Overall F-information')
% hold off
% 
% subplot(2,1,2)
% plot(abs(a'-b'))
% hold on
% plot(b'-0.5)
% plot(abs(b'-0.5)+abs(a'-b'))
%% %% Find the best way to seperate to 4 bands
% close all;
load('transmittance_7_14')
temp_diff=10;
channels=4;
TIME=0.001;
area=5e-6^2;
T_obj = 300; 
T_atm = T_obj-10;
SUBCHANNELS = 1000;
wavelength = linspace(7,14,1000+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,SUBCHANNELS,TIME,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,SUBCHANNELS,TIME,area);
d = 100;
index = round(linspace(1,SUBCHANNELS+1,channels+1));


tem_result_up = zeros(1,1000);
tem_result_down = zeros(1,1000);

for i=1:1000
    c1 = ExpectedPhoton_atm(i);
    c1_record(i) = c1;
    c2 = ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i);
    c2_record(i) = c2;
    [~,trans_index] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index,2);
%         before change to sum Fisher information
%         tem_result(i) = c2^2*transmittance^(2*d)*(log(transmittance))^2/...
%         (c1+c2*transmittance^d); 
    transm_record(i) = transm_single;


    tem_result_up(i) = c2*transm_single^d*log(transm_single);
    tem_result_down(i) = c1+c2*transm_single^d;
end
U = sum(tem_result_up);
D = sum(tem_result_down);

% calculate F-infor of control group
% for ch = 1:channels
%     tem_up = sum(tem_result_up(1,index(ch):index(ch+1)-1))^2;
%     tem_down = sum(tem_result_down(1,index(ch):index(ch+1)-1));
%     control_group(ch) = tem_up/tem_down;
% end

result = zeros(4,1);
index_ijk = zeros(3,1);
count=1;
for i = 2:5:998
    i
    for j= i+1:5:999
        for k = j+1:5:1000
            tem_up = sum(tem_result_up(1,1:i-1))^2;
            tem_down = sum(tem_result_down(1,1:i-1));
            result(1,count) = tem_up/tem_down;

            tem_up = sum(tem_result_up(1,i:j-1))^2;
            tem_down = sum(tem_result_down(1,i:j-1));
            result(2,count) = tem_up/tem_down;

            tem_up = sum(tem_result_up(1,j:k-1))^2;
            tem_down = sum(tem_result_down(1,j:k-1));
            result(3,count) = tem_up/tem_down;
            
            tem_up = sum(tem_result_up(1,k:1000))^2;
            tem_down = sum(tem_result_down(1,j:k-1));
            result(4,count) = tem_up/tem_down;

            index_ijk(:,count) = [i;j;k];

            count = count+1;
        end
    end
end
% tem_result_up = sum(tem_result_up)^2;
% tem_result_down = sum(tem_result_down);
% result = tem_result_up/tem_result_down;

figure()
hold on
plot(sum(result))
[~,biggest_index]=max(sum(result));
% plot([biggest_index,biggest_index],[0.004,0.01])
title('relationship between n and the overall F-information')
xlabel('Position of the interval')
ylabel('Overall F-information')
hold off
%% Find the best way to seperate to multi bands
% close all;
close all;
clear all;
load('transmittance_7_14')
temp_diff=10;
channels=50;
TIME=0.001;
area=5e-6^2;
distance=100;
T_obj = 300; 
T_atm = T_obj-10;
SUBCHANNELS = 1000;
wavelength = linspace(7,14,SUBCHANNELS+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,SUBCHANNELS,TIME,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,SUBCHANNELS,TIME,area);
d = 100;
index = round(linspace(1,SUBCHANNELS+1,channels+1));


tem_result_up = zeros(1,SUBCHANNELS);
tem_result_down = zeros(1,SUBCHANNELS);

for i=1:SUBCHANNELS
    c1 = ExpectedPhoton_atm(i);
    c2 = ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i);
    [~,trans_index] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index,2);
    
    tem_result_up(i) = c2*transm_single^d*log(transm_single);
    tem_result_down(i) = c1+c2*transm_single^d;
end
U = sum(tem_result_up);
D = sum(tem_result_down);

for ch = 1:channels
    tem_up = sum(tem_result_up(1,index(ch):index(ch+1)-1))^2;
    tem_down = sum(tem_result_down(1,index(ch):index(ch+1)-1));
    control_group(ch) = tem_up/tem_down;
end

result = zeros(2,SUBCHANNELS-1);
a = zeros(1,SUBCHANNELS-1);
b = zeros(1,SUBCHANNELS-1);

for ch = 2:30
    ch
%     this index is the default index, starting with [1]
    final_index = 1; 
%     optimize from 2 channels
    for j = 2:(ch)
%         -1 is for debugging
        subresult = zeros(j,1000);
%         the head of the channel will move from the beginning to the end
        for i = 1:1000
            index = final_index;
            if find(index==i)
                continue;
            end
%             find the location of the index just in front of the i
%             insert ind into index
            [~,ind] = max(index(find(index<=i))-i);
            if index(end) <=ind
                index = [index i];
            else
                index = [index(1:ind) i index(ind+1:end)];
            end     
            
            for ch_loop = 1:j
                subresult(ch_loop,i) = 0;
                beg = index(ch_loop);
                if ch_loop==size(index,2)
                    en = 1000;
                else
                    en = index(ch_loop+1);
                end
                tem_up = sum(tem_result_up(1,beg:en))^2;
                tem_down = sum(tem_result_down(1,beg:en));
                subresult(ch_loop,i) = subresult(ch_loop,i) + tem_up/tem_down;
            end
        end
        [~,biggest_i]=max(sum(subresult));
        [~,ind] = max(final_index(find(final_index<=biggest_i))-i);
        if final_index(end) <=ind
            final_index = [final_index biggest_i]
        else
            final_index = [final_index(1:ind) biggest_i final_index(ind+1:end)]
        end  
        
    end
    optimized_result(ch) = sum(FisherInformation(temp_diff, ch,TIME,area,transmittance,final_index));
%     index_1 = round(linspace(1,subchannels,ch+1));
%     referece(ch) = sum(FisherInformation(temp_diff, ch,time,area,transmittance,index_1));
end

% figure()
% plot((optimized_result(find(optimized_result))-referece(find(optimized_result)))./referece(find(optimized_result)))
% title('Improvement of Optimized Channels')
subresult = sum(subresult(:,2:end));
subresult = [0.0045 subresult];
plot(subresult(find(subresult~=0)))

figure()
optimized_result(1) = 0.0045;
plot((optimized_result))
title('Optimized Overall Fisher Information')
%% Using hyperspectal to calculate distance
clear;
run_3_7 = 1;
% #1, we have to prove e^(c2/wavelength/T)>>1, so that we can make
% equation set linear
H = 6.626e-34;
C = 3e8;
K = 1.38e-23;
C2 = 1.438786e-2;
T = 300;
% C2 is the second radiation constant(m*K)
% T is the Kelvin temperature

wavelength = 1e-6:1e-6:20e-6;
exp(C2/T./wavelength)-1;
exp(10e6*H*C./wavelength/K/T)-1;


% #2, Generate a spectral of measurement
load('transmittance_7_14')
global CHANNELS TIME PIXEL_AREA SUBCHANNELS trans_index TEMP_DIFF INTERVAL
TEMP_DIFF=10;
CHANNELS=100;
TIME=0.001;
PIXEL_AREA=5e-6^2;
DISTANCE=600;
T_ATM = 290;
T_OBJ = 290+TEMP_DIFF; 

SUBCHANNELS = 1000;
iter_number = 100;
INTERVAL = 7/SUBCHANNELS;

wavelength = linspace(7,14,SUBCHANNELS+1);
ExpectedPhoton_atm=Expected_Photon(T_ATM,SUBCHANNELS,TIME,PIXEL_AREA);
ExpectedPhoton_obj=Expected_Photon(T_OBJ,SUBCHANNELS,TIME,PIXEL_AREA);

trans_index = zeros(1,CHANNELS);
% Calculate the photon in every subchannel
for i=1:SUBCHANNELS
    [~,trans_index(i)] = min(abs(wavelength(i)-transmittance(:,1)));
    transm_single = transmittance(trans_index(i),2);
    sub_measurement(i) = INTERVAL*(ExpectedPhoton_atm(i)+(ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i))*transm_single^DISTANCE);
end

% Calculate the photon in every channel
global measurement
measurement = zeros(1,CHANNELS);
index = round(linspace(1,SUBCHANNELS+1,CHANNELS+1));
for i=1:CHANNELS
    measurement(i)=sum(sub_measurement(index(i):index(i+1)-1));
end
MEASUREMENT = measurement;

% Plot the measurement
subplot(2,1,1)
bar(sub_measurement)
title('Photons received by subchannels')
subplot(2,1,2)
bar(measurement)
title('Photons received by channels')


% #3 Solve the non-linear system of equations for ideal conditions
fun = @Nonlinear_System;
x0 = [270,280,1];
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
x = fsolve(fun,x0,opts);
disp(['The calculated T_atm is ',int2str(x(1))]);
disp(['The calculated T_obj is ',int2str(x(2))]);
disp(['The calculated distance is ',int2str(x(3))]);
disp('');


% #4 Solve the non-linear system of equations for real conditions
for i=1:CHANNELS
    measurement(i) = poissrnd(MEASUREMENT(i));
end
fun = @Nonlinear_System;
x0 = [270,270,0];
opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
x = fsolve(fun,x0,opts);

disp(['The calculated T_atm is ',int2str(x(1))]);
disp(['The calculated T_obj is ',int2str(x(2))]);
disp(['The calculated distance is ',int2str(x(3))]);
disp('');



% #5 Test the results when know T_atm
x0 = [270,0];
for test_itr = 1:iter_number
    disp(['iteration: ',int2str(test_itr)]);
    for i=1:CHANNELS
    measurement(i) = poissrnd(MEASUREMENT(i));
    end
    x = fsolve(fun,x0,opts);
    test_T(test_itr) = x(1);
    test_dist(test_itr) = x(2);
end

figure()
subplot(2,1,1)
hist(test_T,12)
title('Calculated object temperature(T=300K)')

subplot(2,1,2)
hist(test_dist,12)
title('Calculated object distance(d=600)')


%  #6 Test the results when T_atm unknown
x0 = [270,270,0];
for test_itr = 1:iter_number
    disp(['iteration: ',int2str(test_itr)]);
    for i=1:CHANNELS
        measurement(i) = poissrnd(MEASUREMENT(i));
    end
    
    x = fsolve(fun,x0,opts);
    test_T_atm(test_itr) = x(1);
    test_T_obj(test_itr) = x(2);
    test_dist(test_itr) = x(3);
end
figure()
subplot(3,1,1)
hist(test_T_atm,12)
title('Calculated atmosphere temperature(T=290K)')

subplot(3,1,2)
hist(test_T_obj,12)
title('Calculated object temperature(T=300K)')

subplot(3,1,3)
hist(test_dist,12)
title('Calculated object distance(d=600)')

%  #7 Plot the Fisher information of every channel
figure()
F_inform = FisherInformation(TEMP_DIFF, CHANNELS,DISTANCE,TIME,PIXEL_AREA,transmittance,index);
bar(F_inform)
title('Fisher information distribution')

% %  Plot the Fisher information of every channel in different distance
% 2m~200m
% F_inform = zeros(100,100);
% for i=1:100
%     disp(['iteration: ',int2str(i)]);
%     DISTANCE = i*2;
%     F_inform(i,:) = FisherInformation(TEMP_DIFF, CHANNELS,DISTANCE,TIME,PIXEL_AREA,transmittance,index);
% end
% 
% figure()
% bar3(F_inform)
% title('Fisher information distribution varies with distance')



%  #8 Use optimization of negative log instead of solve nonlinear system
cost_fun = @Negative_Log;
x0 = [270,280,200];
% A = [-1,0,0;0,-1,0;0,0,-1];
% b = [0,0,0];
A = [];
b = [];
for test_itr = 1:iter_number
    disp(['iteration: ',int2str(test_itr)]);
    for i=1:CHANNELS
    measurement(i) = poissrnd(MEASUREMENT(i));
    end
    optimized_x = fmincon(cost_fun,x0,A,b);
    Tatm_unKnown_Tatm(test_itr) = optimized_x(1);
    Tobj_unKnown_Tatm(test_itr) = optimized_x(2);
    Distance_unKnown_Tatm(test_itr) = optimized_x(3);
    Likelihood_fun(test_itr) = Negative_Log(optimized_x);
    Likelihood(test_itr) = Negative_Log([290,300,200]);
end

figure()
subplot(3,1,1)
hist(Tatm_unKnown_Tatm,12)
title('Calculated atm temperature(T=290K) with Neglog')

subplot(3,1,2)
hist(Tobj_unKnown_Tatm,12)
title('Calculated object temperature(T=300K) with Neglog')

subplot(3,1,3)
hist(Distance_unKnown_Tatm,12)
title(['Calculated object distance(d=',int2str(DISTANCE),'m) with Neglog'])


%% Explore CRB of different distance, temperature and time
clear distance CRB_n30 CRB_n20 CRB_n10 CRB_p10 CRB_p20 CRB_p30 CRB_p100;
distance = 100:100:2000;

% CRB for different distance and temperature
t_obj = T_ATM+[-30,-20,-10,10,20,30,40];

for i =1:size(distance,2)
    disp(['CRB_T iteration: ',int2str(i)]);
    DISTANCE = distance(i);
    CRB_n30(i) = CramerRaoBound(T_ATM,t_obj(1),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_n20(i) = CramerRaoBound(T_ATM,t_obj(2),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_n10(i) = CramerRaoBound(T_ATM,t_obj(3),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_p10(i) = CramerRaoBound(T_ATM,t_obj(4),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_p20(i) = CramerRaoBound(T_ATM,t_obj(5),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_p30(i) = CramerRaoBound(T_ATM,t_obj(6),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
    CRB_p40(i) = CramerRaoBound(T_ATM,t_obj(7),DISTANCE,CHANNELS,SUBCHANNELS,TIME,PIXEL_AREA);
end

figure()
hold on
plot(distance,CRB_n30,'-o','Color',[0,0.7,0.9])
plot(distance,CRB_n20,'-o')
plot(distance,CRB_n10,'-o')
plot(distance,CRB_p10,'-o')
plot(distance,CRB_p20,'-o')
plot(distance,CRB_p30,'-o')
plot(distance,CRB_p40,'-o')
hold off

legend('260K','270K','280K','300K','310K','320K','330K')
title('Cramer-Rao bound for a blackbody at 290K atmosphere')
xlabel('distance')
ylabel('CRLB')


% CRB for different distance and temperature
time = [0.001 0.005 0.01 0.05 0.1];
T_OBJ = T_ATM+10;

for i =1:size(distance,2)
    disp(['CRB_TIME iteration: ',int2str(i)]);
    DISTANCE = distance(i);
    CRB_1ms(i) = CramerRaoBound(T_ATM,T_OBJ,DISTANCE,CHANNELS,SUBCHANNELS,time(1),PIXEL_AREA);
    CRB_5ms(i) = CramerRaoBound(T_ATM,T_OBJ,DISTANCE,CHANNELS,SUBCHANNELS,time(2),PIXEL_AREA);
    CRB_10ms(i) = CramerRaoBound(T_ATM,T_OBJ,DISTANCE,CHANNELS,SUBCHANNELS,time(3),PIXEL_AREA);
    CRB_50ms(i) = CramerRaoBound(T_ATM,T_OBJ,DISTANCE,CHANNELS,SUBCHANNELS,time(4),PIXEL_AREA);
    CRB_100ms(i) = CramerRaoBound(T_ATM,T_OBJ,DISTANCE,CHANNELS,SUBCHANNELS,time(5),PIXEL_AREA);
end

figure()
hold on
plot(distance,CRB_1ms,'-o')
plot(distance,CRB_5ms,'-o')
plot(distance,CRB_10ms,'-o')
plot(distance,CRB_50ms,'-o')
plot(distance,CRB_100ms,'-o')
hold off

legend('1ms','5ms','10ms','50ms','100ms')
title('Cramer-Rao bound for a 300K blackbody at 290K atmosphere')
xlabel('distance')
ylabel('CRLB')
%% Generate the non-linear system of equations
function F = Nonlinear_System(x)
global measurement CHANNELS TIME PIXEL_AREA SUBCHANNELS trans_index INTERVAL
if size(x,2)==2
    T_ATM = 290;
    T_OBJ = x(1);
    DISTANCE = x(2);
else
    T_ATM = x(1);
    T_OBJ = x(2);
    DISTANCE = x(3);
end
F = zeros(1,CHANNELS);
load('transmittance_7_14');


ExpectedPhoton_atm=Expected_Photon(T_ATM,SUBCHANNELS,TIME,PIXEL_AREA);
ExpectedPhoton_obj=Expected_Photon(T_OBJ,SUBCHANNELS,TIME,PIXEL_AREA);
index = round(linspace(1,SUBCHANNELS+1,CHANNELS+1));
for i=1:SUBCHANNELS
    transm_single = transmittance(trans_index(i),2);
    sub_measurement(i) = INTERVAL*(ExpectedPhoton_atm(i)+(ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i))*transm_single^DISTANCE);
end
for i=1:CHANNELS
    F(i)=sum(sub_measurement(index(i):index(i+1)-1))-measurement(i);
end
end

function Log_Likelihood = Negative_Log(x)
global measurement CHANNELS TIME PIXEL_AREA SUBCHANNELS trans_index TEMP_DIFF INTERVAL
if size(x,2)==2
    T_ATM = 300-TEMP_DIFF;
    T_OBJ = x(1);
    DISTANCE = x(2);
else
    T_ATM = x(1);
    T_OBJ = x(2);
    DISTANCE = x(3);
end
eta = zeros(1,CHANNELS);
load('transmittance_7_14');


ExpectedPhoton_atm=Expected_Photon(T_ATM,SUBCHANNELS,TIME,PIXEL_AREA);
ExpectedPhoton_obj=Expected_Photon(T_OBJ,SUBCHANNELS,TIME,PIXEL_AREA);
index = round(linspace(1,SUBCHANNELS+1,CHANNELS+1));
for i=1:SUBCHANNELS
    transm_single = transmittance(trans_index(i),2);
    sub_measurement(i) = INTERVAL*(ExpectedPhoton_atm(i)+(ExpectedPhoton_obj(i)-ExpectedPhoton_atm(i))*transm_single^DISTANCE);
end
for i=1:CHANNELS
    eta(i)=sum(sub_measurement(index(i):index(i+1)-1));
end
Log_Likelihood = sum(eta-measurement.*log(eta));
end





