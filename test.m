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
time=0.001;
area=5e-6^2;
distance=100;
expt = Expected_Photon(temp_diff,channels,time,area);

% Fisher information of 10 channels
wavelength = linspace(7,14,channels);
figure(1)
subplot(2,2,1)
result_10 = FisherInformation(temp_diff, 10,time,...
    area,transmittance);
plot(result_10)
title('Fisher information of 10 channels')
% Fisher information of 100 channels
subplot(2,2,2)
result_100 = FisherInformation(temp_diff, 100,time,...
    area,transmittance);
plot(result_100)
title('Fisher information of 100 channels')
% Fisher information of 500 channels
subplot(2,2,3)
result_500 = FisherInformation(temp_diff, 500,time,...
    area,transmittance);
plot(result_500)
title('Fisher information of 500 channels')

subplot(2,2,4)
plot(transmittance(:,1),1-transmittance(:,2))
title('1-transmittance')

% study the relationship between channels and total F-information
figure(2)
num_channel=1:100;
for i=1:100
    result(i) = sum(FisherInformation(temp_diff,num_channel(i),time,...
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
    result(i) = sum(FisherInformation(diff(i),num_channel,time,...
    area,transmittance));
end
plot(diff,result)
title('Relationship between temperature difference and F-information')


%% Find the best way to seperate the band
close all;
load('transmittance_7_14')
temp_diff=10;
channels=2;
time=0.001;
area=5e-6^2;
distance=100;
T_obj = 300; 
T_atm = T_obj-10;
wavelength = linspace(7,14,3000+1);
ExpectedPhoton_atm=Expected_Photon(T_atm,3000,time,area);
ExpectedPhoton_obj=Expected_Photon(T_obj,3000,time,area);
d = 100;
index = round(linspace(1,3001,channels+1));


tem_result_up = zeros(1,3000);
tem_result_down = zeros(1,3000);

for i=1:3000
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
result = zeros(2,2999);
a = zeros(1,2999);
b = zeros(1,2999);
for i = 1:2999
    tem_up = sum(tem_result_up(1,1:i))^2;
    tem_down = sum(tem_result_down(1,1:i));
    result(1,i) = tem_up/tem_down;
    a(1,i) = tem_down;
    b(1,i) = sum(tem_result_up(1,1:i));
    
    tem_up = sum(tem_result_up(1,i+1:3000))^2;
    tem_down = sum(tem_result_down(1,i+1:3000));
    result(2,i) = tem_up/tem_down;
    ratio_a(i) = D*b(i)/U;
    ratio_a2(i) = D*b(i)/(2*b(i)-U);
    ratio_b(i) = a(i)*U/D;
end
% tem_result_up = sum(tem_result_up)^2;
% tem_result_down = sum(tem_result_down);
% result = tem_result_up/tem_result_down;

figure(1)
subplot(3,1,1)
hold on
plot(sum(result))
plot([1,2999],[0.0046,0.0046],'r')
[~,biggest_index]=max(sum(result));
plot([biggest_index,biggest_index],[0.004,0.01])
hold off

subplot(3,1,2)
[~,biggest_diff_index]=max(abs(ratio_a-a));
hold on
plot(ratio_a')
plot(a')
plot([biggest_diff_index,biggest_diff_index],[0,6e5])
plot(ratio_a2')
legend('ratio a')
legend('a')
hold off

subplot(3,1,3)
[~,biggest_diff_index]=max(abs(ratio_b-b));
plot(ratio_b')
legend('ratio b')
hold on
plot(b')
legend('b')
plot([biggest_diff_index,biggest_diff_index],[0,-50])
hold off

