function [ Fingerprint_R, Fingerprint_G, Fingerprint_B ] = GetRGBFingerprint( fileName, fileName2, fileName3 )

Images(1).name = fileName;  
if nargin == 2
    Images(2).name = fileName2;  
elseif nargin == 3
    Images(2).name = fileName2;  
    Images(3).name = fileName3; 
end

RP = getFingerprint(Images);
Fingerprint_R = WienerInDFT(RP(:,:,1),std2(RP(:,:,1)));
Fingerprint_G = WienerInDFT(RP(:,:,2),std2(RP(:,:,2)));
Fingerprint_B = WienerInDFT(RP(:,:,3),std2(RP(:,:,3)));

end

