function [ suspicious_rgb_imx ] = Modify_hue_degree( image, degree )

suspicious_hsv_imx = rgb2hsv(image);
suspicious_hue_imx = suspicious_hsv_imx(:,:,1)*360; % extract hue
suspicious_hue_imx = mod((suspicious_hue_imx + degree), 360);
suspicious_hsv_imx(:,:,1) = suspicious_hue_imx/360;
suspicious_rgb_imx = hsv2rgb(suspicious_hsv_imx);
%imwrite(uint8(suspicious_rgb_imx),savefn);

end

