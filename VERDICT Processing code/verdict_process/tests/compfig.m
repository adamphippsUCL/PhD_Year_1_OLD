function compfig(Aref, B, im_rng, diff_rng )
% COMPFIG Compare figures. Used in testing
%
% compfig(Aref, B, im_rng, diff_rng )

dbs = dbstack ;

if length(dbs) > 1
    str = dbs(2).name ;
else
    str = '';
end

hf = figure('Name',['compfig: ', str]) ;
t = tiledlayout('flow','Parent',hf,'Padding','tight','TileSpacing','tight');
ax1 = nexttile(t) ;
imshow(Aref,im_rng)
xlabel('ref')
ax2 = nexttile(t) ;
imshow(B,im_rng)
ax3 = nexttile(t) ;
imdiff = mat2gray(B-Aref, diff_rng) ;
[imdiffind, cmap] = gray2ind(imdiff) ;
cmap(1,:) = [1 0 0] ;
cmap(end,:) = [1 0 0] ;
imRGB = ind2rgb(imdiffind, cmap) ;
imshow(imRGB)
xlabel(['range: [',num2str(diff_rng(1)), ' ', num2str(diff_rng(2)),']'])
linkaxes([ax1 ax2 ax3])

end