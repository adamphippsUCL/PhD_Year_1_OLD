function dmap = explore_vcmap(fIC,v)
% EXPLORE_VCMAP Explore colormap for VERDICT processed data (fIC maps)
%
% Would like to convey, anatomy, fIC (low and high) and uncertainty.
% Diverging colormap with white point adjusted to Singh et al 0.42 looks OK.
% though no anatomy.
%
% Use of product (not transpency) blending inpsired by geography 
% (relief shading and isoluminace maps), but tends not to look great here, 
% even with some manipulation of
% base image. I think the problem is the variety of gray levels.
%
% Previous use of alpha = fIC, fixed colour is not bad (but no certainty)
%
%
% See https://colorcet.com and https://doi.org/10.48550/arXiv.1509.03700
% Peter Kovesi. Good Colour Maps: How to Design Them.
% arXiv:1509.03700 [cs.GR] 2015
%
%

testFile = 'colourmaptest.tif' ;

[A, map] = imread(testFile) ;

% D11 is isoluminant, R3 is divergent rainbow, D09 divergent (also for
% multiplication)
d01map=colorcet('D09') ;

% Skew colormap to put mid-white at fic=0.42 for divergent
whiteVal = 0.42 ;

cmv = linspace(0,1,size(d01map,1)) ;
cmq = cmv/whiteVal*0.5 ;

extrapVal = d01map(end,:) ;

redv   = interp1(cmv, d01map(:,1), cmq,"linear",extrapVal(1)) ;
greenv = interp1(cmv, d01map(:,2), cmq,"linear",extrapVal(2)) ;
bluev  = interp1(cmv, d01map(:,3), cmq,"linear",extrapVal(3)) ;

dmap = cat(2, redv(:), greenv(:), bluev(:)) ;

figure(Name='test image')
tiledlayout(2,1)
nexttile
imshow(A, dmap)

nexttile
lw=2;
plot(dmap(:,1),'r','LineWidth',lw), hold on
plot(dmap(:,2),'g','LineWidth',lw)
plot(dmap(:,3),'b','LineWidth',lw)
axis([0 256 0 1])
grid


%
v = mat2gray(v,[0 prctile(v(:),98)]) ;

fIC(fIC<0)=0;
fIC(fIC>1)=1;
nc = size(dmap,1) ;

[ny,nx,nz] = size(fIC) ;
ficovRGB = zeros([ny nx nz 3]) ;
ficRGB = zeros([ny nx nz 3]) ;

for iz = 1:nz
    indFIC = round(double(fIC(:,:,iz))*(nc-1) + 1) ; % fIC is [0 1], convert to [1 nc]
    overRGB = ind2rgb(indFIC,dmap) ;
    fsl = imadjust(v(:,:,iz),[],[0.3 1],1) ; % to prevent regions that are too dark
    ficRGB(:,:,iz,:) = reshape(overRGB,[ny nx 1 3]);

    ficovRGB(:,:,iz,:) = cat(4, fsl.*overRGB(:,:,1), fsl.*overRGB(:,:,2), fsl.*overRGB(:,:,3) ) ;
end

sviewer(ficovRGB,[],isrgb=true)
sviewer(ficRGB,[],isrgb=true)
