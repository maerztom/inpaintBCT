% images
im = double(imread('newOrleans.png'));
immask = double(imread('newOrleans_mask.png'));

% mask
mask = getMask(im,immask);

% inpaint
tic;
imr = inpaintBCT(immask,'orderD',mask,'guidanceC',[4 25 2 3]);
toc

% show result
figure; 
rect = get(gcf,'OuterPosition');
rect(1) = rect(1)-round(rect(3)/2);
rect(3) = 2*rect(3);
set(gcf,'OuterPosition',rect);
clear rect
subplot(1,2,1), subimage(uint8(im)), axis off, axis image, title('original')
subplot(1,2,2), subimage(uint8(imr)), axis off, axis image, title('inpainted')