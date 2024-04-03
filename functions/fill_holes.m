function [chl2] = fill_holes(chl,size,inp_meth)
% Luke Carberry, 03/16/2022
% function to interpolate over small holes
% optimal hole size is 14
% requires inpaint_nans from John D'Errico 
% 1: chl = Landsat chlorophyll matrix
% 2: size = upper size cutoff (in pixels) for small holes. i.e. the smallest
% hole size that escapes interpolation, so size = 1 returns the same matrix
% 3: inpaint_nans method, 5 is nearest neighbor, 1 is least squares. Visual
% inspect for which one looks most natural. My preference is 5, due to
% simplicity

%uncomment the imshow commands to see what each step does

chlmask = chl;
I = ~isnan(chl);
chlmask(I) = 1;
J = isnan(chl);
chlmask(J) = 0;
% imshow(chlmask')

chl_fill = imfill(chlmask,8);
% imshow(chl_fill')
% title('All holes filled')

holes = chl_fill & ~chlmask;
% imshow(holes')
% title('Hole pixels identified')

bigholes = bwareaopen(holes,size,4);
% imshow(bigholes')
% title('Only the big holes')

smallholes = holes & ~bigholes;
% imshow(smallholes')
% title('Only the small holes')

% new = chlmask | smallholes;
% imshow(new')
% title('Small holes filled')

chl_inp = chl;
chl_inp(isnan(chl)) = mean(chl,'all','omitnan');
chl_inp(smallholes) = NaN;
chl_inp = inpaint_nans(chl_inp,inp_meth);
chl2 = chl;
chl2(smallholes) = chl_inp(smallholes);

% imshow(chl2')
% title('Small holes filled')
