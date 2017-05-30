clear all
close all 
clc

disp('1 step : Extract R,G,B reference patterns');
fn_rf_1 = 'Images\P1.jpg';
fn_rf_2 = 'Images\P2.jpg';
fn_rf_3 = 'Images\P3.jpg';
[ Fingerprint_R, Fingerprint_G, Fingerprint_B ] = Get_RGB_Fingerprint(fn_rf_1, fn_rf_2, fn_rf_3);

disp('2 step : Make suspicious Image');
fn_ori = 'Images\P1.jpg';
original_img = double(imread(fn_ori));
degree = 160;
suspicious_rgb_imx = Modify_hue_degree( original_img, degree );

disp('3. Make hue shifted images from 0 to 359 degree.');
disp('4. Extract R, G, B noise residual of shifted images.');
disp('5. Calculate correlation'); 
detection = zeros([360,1]);
for i = (0:1:359)
    % 3. Make hue shifted images from 0 to 359 degree.
    hue_shifted_img = Modify_hue_degree( suspicious_rgb_imx, i);
    % 4. Extract R, G, B noise residual of shifted images.
    [Noise_R, Noise_G, Noise_B] = Extract_RGB_Noise(hue_shifted_img);
    % 5. Calculate correlation 
    C = crosscorr(Noise_R,Fingerprint_R) + crosscorr(Noise_G,Fingerprint_G) + crosscorr(Noise_B,Fingerprint_B);
    detection(i+1) = getfield(PCE(C),'PCE');
end

disp('6. Find degree.');
[max_corr,degree] = max(detection);
disp('Shifted degree :');
degree = mod(360 - (degree - 1) ,360);
disp(degree);

disp('7. Restored image.');
restored_rgb_imx = Modify_hue_degree( suspicious_rgb_imx, - degree );

figure('units','pixels','position',[200 200 1200 400]);
axes('units','norm','outerposition',[0 0 0.33 1],'position',[0 0 0.33 1]) 
imshow(uint8(suspicious_rgb_imx)); title('suspicious image');
axes('units','norm','outerposition',[0.33 0 0.33 1],'position',[0.33 0 0.33 1])
imshow(uint8(restored_rgb_imx)); title('restored image');
axes('units','norm','outerposition',[0.66 0 0.33 1],'position',[0.66 0 0.33 1])
imshow(uint8(original_img)); title('original image');





