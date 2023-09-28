% MATLAB script to display gridsearch results

% Displaying dSdR over b value plane with highest Sensitivity
fig1 = figure;
pcolor(Ds(:,:,I3), ds(:,:,I3), Sensitivity_array(:,:,I3));
colormap("jet");
colorbar;
xlabel('Delta');
ylabel('delta');
title(['Sensitivity at R = ' num2str(Rcell) ', b = ' num2str(bs(I3)) ' ADC = ' num2str(ADC)]);
saveas(fig1, ['Figures/grid search/Sensitivity at R = ' num2str(Rcell) ', b = ' num2str(bs(I3)) ' ADC = ' num2str(ADC) '.png'] )

% Displaying signal over b value plane with highest Sensitivity
fig2= figure;
pcolor(Ds(:,:,I3), ds(:,:,I3), Signal_array(:,:,I3));
colormap("jet");
colorbar;
xlabel('Delta');
ylabel('delta');
title(['Signal at R = ' num2str(Rcell) ', b = ' num2str(bs(I3)) ' ADC = ' num2str(ADC)]);
saveas(fig2, ['Figures/grid search/Signal at R = ' num2str(Rcell) ', b = ' num2str(bs(I3)) ' ADC = ' num2str(ADC) '.png'])

% Displaying dSdR over delta value plane with highest Sensitivity
fig3 = figure;
pcolor(squeeze( Ds(:,I2,:) ), squeeze( Bs(:,I2,:) ), squeeze( Sensitivity_array(:,I2,:)));
colormap("jet");
colorbar;
xlabel('Delta');
ylabel('b');
title(['Sensitivity @ R = ' num2str(Rcell) ', delta = ' num2str(deltas(I2)) ' ADC = ' num2str(ADC)]);
saveas(fig3, ['Figures/grid search/Sensitivity at R = ' num2str(Rcell) ', delta = ' num2str(deltas(I2)) ' ADC = ' num2str(ADC) '.png'])

% Displaying signal over delta value plane with highest Sensitivity
fig4 = figure;
pcolor(squeeze( Ds(:,I2,:) ), squeeze( Bs(:,I2,:) ), squeeze( Signal_array(:,I2,:)));
colormap("jet");
colorbar;
xlabel('Delta');
ylabel('b');
title(['Signal at R = ' num2str(Rcell) ', delta = ' num2str(deltas(I2)) ' ADC = ' num2str(ADC) ]);
saveas(fig4, ['Figures/grid search/Signal at R = ' num2str(Rcell) ', delta = ' num2str(deltas(I2)) ' ADC = ' num2str(ADC) '.png'])

% Displaying dSdR over Delta value plane with highest Sensitivity
fig5 = figure;
pcolor(squeeze( ds(I1,:,:) ), squeeze( Bs(I1,:,:) ), squeeze( Sensitivity_array(I1,:,:)));
colormap("jet");
colorbar;
xlabel('delta');
ylabel('b');
title(['Sensitivity at R = ' num2str(Rcell) ', Delta = ' num2str(Deltas(I1)) ' ADC = ' num2str(ADC)]);
saveas(fig5, ['Figures/grid search/Sensitivity at R = ' num2str(Rcell) ', Delta = ' num2str(Deltas(I1)) ' ADC = ' num2str(ADC) '.png'] )

% Displaying signal over Delta value plane with highest Sensitivity
fig6 = figure;
pcolor(squeeze( ds(I1,:,:) ), squeeze( Bs(I1,:,:) ), squeeze( Signal_array(I1,:,:)));
colormap("jet");
colorbar;
xlabel('delta');
ylabel('b');
title(['Signal at R = ' num2str(Rcell) ', Delta = ' num2str(Deltas(I1)) ' ADC = ' num2str(ADC)]);
saveas(fig6, ['Figures/grid search/Signal at R = ' num2str(Rcell) ', Delta = ' num2str(Deltas(I1)) ' ADC = ' num2str(ADC) '.png'])