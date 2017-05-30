function [ Noise_R, Noise_G, Noise_B ] = Extract_RGB_Noise ( image )

    Noisex = NoiseExtractFromImage(image,3,'color');
    Noise_R = WienerInDFT(Noisex(:,:,1),std2(Noisex(:,:,1)));
    Noise_G = WienerInDFT(Noisex(:,:,2),std2(Noisex(:,:,2)));
    Noise_B = WienerInDFT(Noisex(:,:,3),std2(Noisex(:,:,3)));

end

