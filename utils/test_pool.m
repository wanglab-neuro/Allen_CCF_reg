% imref = TVL1denoise(b, 0.1);
% imcur = TVL1denoise(a, 0.1);
% imref = TVL1denoise(normalize(b .* bw),0.1);
% imcur = TVL1denoise(a,0.1);
% imref = normalize(b .* bw);
% imcur = a;
% imref = imgaussfilt(normalize(b .* bw), 3);
% imcur = imgaussfilt(a, 3);
% imref = imgaussfilt(normalize(b), 1);
imref = imgaussfilt(imgref, 9);
imcur = imgaussfilt(imgcur, 9);
imcur(imcur < 0.1) = 0;
imref(imref < 0.1) = 0;

[~,img,~]=klt_ref_track(imcur,imref);

imref = TVL1denoise(imref, 0.1);
imcur = TVL1denoise(img, 0.1);

imref = double(imref > 0.1);
imcur = double(imcur > 0.1);

[imgt, sx, sy] = logdemons_unit(imcur, imref, [], [], 2, 2, 1);
[D,img]=imregdemons(imcur,imref,200,'AccumulatedFieldSmoothing',2,'displaywaitbar',false);

imgg = normalize(double(imresize(aa(:, :, 2), size(imref))));
imgg = iminterpolate(imgg, D(:, :, 1), D(:, :, 2));


[~,tps,sps]=klt2(imref,imcur,[],0.01);
tform = fitgeotrans(tps,sps,'lwm',40);
img = imwarp(imcur,tform,'smoothedges',true);
imagesc(img)
showMatchedFeatures(imref,imcur,sps,tps)


imgg = normalize(double(imresize(aa(:, :, 2), size(imref))));
for i = 1: length(sx)
    imgg = iminterpolate(imgg, sx{i}, sy{i});
end

[optimizer,metric] = imregconfig('multimodal');
imgt = imregister(imcur,imref,'affine',optimizer,metric);

FileName = 'E:\Jinghao\for_others\Jun\allenCCF\P56_atlasVolume\atlasVolume.raw';
% 25 micron volume size
size = [528 320 456];
% VOL = 3-D matrix of atlas Nissl volume
fid = fopen(FileName, 'r', 'l' );
VOL = fread( fid, prod(size), 'uint8' );
fclose( fid );
VOL = reshape(VOL,size);










