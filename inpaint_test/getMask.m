function mask = getMask(im1,im2)
mask = [];

ND1 = ndims(im1);
if (ND1 > 3) | ( ND1 < 2)
    disp('First variable is not an image');
    return;
end
ND2 = ndims(im2);
if (ND2 > 3) | ( ND2 < 2)
    disp('Second variable is not an image');
    return;
end
 
[N1,M1,C1] = size(im1);
[N2,M2,C2] = size(im2);

if (N1 ~= N2) | (M1 ~= M2) | (C1 ~= C2)
    disp('Dimension mismatch');
    return;
end

mask = zeros(N1,M1);
diffim = im1-im2;

for k=1:C1
    ind = find(diffim(:,:,k));
    mask(ind) = 1;
end


