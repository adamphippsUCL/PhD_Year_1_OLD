function T1display(T1)
% T1DISPLAY  Display T1 values
% T1display(T1)
% T1 can be 2D or 3D.
%
% Copyright 2020-2021. David Atkinson  University College London
% D.Atkinson@ucl.ac.uk
%
% See also IRT1

if size(T1,3) > 1
    is3D = true;
else
    is3D = false;
end

prompt={'Enter minimum T1:', 'Enter max T1:', 'Figure Name/Title:'};
name='Input T1 display';
numlines=1;
defaultanswer={'200','3000',inputname(1)};
answer = inputdlg(prompt,name,numlines,defaultanswer) ;
lb = str2num(answer{1}) ;
ub = str2num(answer{2}) ;
figname = answer{3} ;

T1(isnan(T1)) = 0 ;

hf = figure('Name',figname) ;
if verLessThan('matlab','8.4.0')
    colormap jet
else
    colormap parula
end
cmap = colormap ;

if is3D
    T1 = reshape(T1,[size(T1,1) size(T1,2) 1 size(T1,3)]) ;
    montage(T1,'DisplayRange',[lb ub])
else
    
    imshow(T1,[lb ub],'Border','loose','Colormap',cmap)
end


colorbar
title(figname)
disp(['lb ',num2str(lb),' ub: ',num2str(ub),' name: ', ...
    figname,' Opened figure: ',num2str(fignum(hf))])

