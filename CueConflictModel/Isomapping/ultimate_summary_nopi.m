function ultimate_summary_nopi(param_two,param_one, vis2_loc, total_cells, polar, summaries)
%N.B: timsteps = simulation time (in seconds) x 100

correlations = zeros(1,numel(param_two));
pSD = zeros(1,numel(param_two));
discrepancy = zeros(numel(param_two), numel(vis2_loc)+1);
discmean = zeros(numel(param_two), numel(vis2_loc)+1);
discunder = zeros(numel(param_two), numel(vis2_loc)+1);

for adx = 1:numel(param_two)
    dirstringW = [num2str(param_two(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/_no_PI/parameter_search/MultiThread'];
    tier_1_path = [parentpath,'/deg',dirstringW];    %Change depending on what parameter we're plotting.
 
    cd(tier_1_path); 
    
    for idx = 1:numel(param_one)
        dirstringS = [num2str(param_one(idx))];
        tier_2_path = [tier_1_path,'/',dirstringS];
    
        cd(tier_2_path); 
        
        Pvectors = zeros(1,numel(vis2_loc)+1);
        Pvectors(1) = 0;

        undershoot = zeros(1,numel(vis2_loc)+1);
        undershoot(1) = 0;
        
        Means = zeros(1,numel(vis2_loc));
        ResVecLens = zeros(1,numel(vis2_loc));
        Vars = zeros(1,numel(vis2_loc));
        SDs = zeros(1,numel(vis2_loc));
        
        
         for jdx = 1:numel(vis2_loc)
    
            dirstringL = [num2str(vis2_loc(jdx))];
            tier_3_path = [tier_2_path, '/', dirstringL];
            cd(tier_3_path);
            
            fid = fopen('Pvector.dat', 'r');
            str = fgets(fid);
            Pvectors(jdx+1) = sscanf(str, 'Conflict population vector is: %f');

            undershoot(jdx+1) = vis2_loc(jdx) - Pvectors(jdx+1);
            discunder(adx, jdx+1) = vis2_loc(jdx) - Pvectors(jdx+1);
            Pvectors(jdx+1) = Pvectors(jdx+1) - 180;

            
            fclose(fid);           
            
           
            
            final_vector = zeros(1,total_cells);

vector_compare = zeros(1,total_cells);

%Get final ff_weight vectors

fid = fopen('final_ff_weight_vectors.bdat', 'rb');

final_vector = fread(fid, 'float32');

fclose(fid);


for cell = 1:total_cells

    cell_dir=(cell*0.72);

    vector_compare(cell) = (final_vector(cell)- cell_dir);

    %if vector_compare(cell) < 0
     %   vector_compare(cell) = abs(vector_compare(cell));
    %end

    if abs(vector_compare(cell)) > 180
        vector_compare(cell) = 360 - abs(vector_compare(cell));
    end
   
end

sine_sum = 0;
cos_sum = 0;

for cell = 1: numel(vector_compare)
    
    sine_sum=sine_sum+sind(vector_compare(cell));
    cos_sum=cos_sum+cosd(vector_compare(cell));
    
end

if(sine_sum>0&&cos_sum>0)
    
    mean_cells=atand(sine_sum/cos_sum);
    
elseif(cos_sum<0)
    
    mean_cells=atand(sine_sum/cos_sum)+180;
    
else
    
    mean_cells=atand(sine_sum/cos_sum)+360;
    
end



%NOW CONVERT TO RADIANS TO WORK OUT CIRCULAR VARIANCE AND S.D.


rad_veccomp = zeros(1, total_cells);
rad_veccomp(:) = (abs(vector_compare)) * (pi/180); %This is the array of circular data points I am working out the variance of, converted to radians.

%Plotting a rose plot for the data

if polar
bin = 0:10:total_cells;
rose2(rad_veccomp, bin);
title('HD Layer Weight Shift','FontSize', 24);
set(gcf,'Position', get(0,'Screensize'));   %Maximise figure to look good when saved.
set(gcf, 'PaperPositionMode', 'auto');      %Overwite tendency of 'saveas' command to resize figure back again.
saveas(gcf, 'RosePlot', 'epsc');
close(gcf);
end

%Doing part of transformation of vector_compare into 2d vectors
a = sum(sin(rad_veccomp(:)));
b = sum(cos(rad_veccomp(:)));


%Calculating mean resultant vector
mean_res_vec = zeros(2,1);
mean_res_vec(1,1) = (a/total_cells); 
mean_res_vec(2,1) = (b/total_cells); %N.B. Divided by number of cells (number of data points) as per vector averaging to get mean res vec.

%Calculating resultant vector length
res_vec_len = norm(mean_res_vec); %N.B. According to Berens (2009), this resultant vector length measures spread, closer to 1 is more concentrated on mean direction.

%Now calculating mean, variance and sd
poss_mean = atan2((a/total_cells),(b/total_cells)); %This is the same as the mean calculated above, was just checking.
variance = 1-res_vec_len;
sd = sqrt(2* variance);
poss_mean = abs(poss_mean);
         
 
newfile = fopen('cell_shifts.txt', 'w');
for idx = 1: total_cells
    fprintf(newfile, '%f\n', abs(vector_compare(idx)));
end
fclose(newfile);

fileID = fopen('_standdev.txt','w');
fprintf(fileID, 'Mean = %f \n Resultant Vector Length = %f\n Variance = %f\n S.D = %f\n',(poss_mean*(180/pi)), res_vec_len, variance, (sd * (180/pi)));
fclose(fileID);

ResVecLens(jdx) = res_vec_len;
Vars(jdx) = variance;
SDs(jdx) = sd;
Means(jdx) = poss_mean;

discmean(adx, jdx+1) = (poss_mean * (180/pi));            
            
            
            cd(tier_2_path);
            
         end
         
         
         
         
         x = [0:10:180];
         
         
         
         figure();
         plot(x,undershoot, 'b', 'LineWidth', 2);
         set(gca, 'FontSize',24);
         %hold on;
         %plot([0,180],[0,180],'k', 'LineWidth', 2);
         hold on
         plot(x,discmean(adx,:), '--r', 'LineWidth', 2);
         
         xlabel('Conflict (Deg)', 'FontSize', 32);
         ylabel('Undershoot/Remap(Deg)', 'FontSize', 32);
         
         ylim([0,180]);
         xlim([0,180]);
         set(gca, 'FontSize',24);
         set(gca, 'YTick',[0:60:180]);
         set(gca, 'XTick',[0:60:180]);
         saveas(gcf,'summary_graph', 'epsc');
         close(gcf);

         newfile = fopen('PV_summary.dat', 'wt');
         for jdx = 1:(numel(vis2_loc)+1)
             fprintf(newfile, '%f', Pvectors(jdx));
             fprintf(newfile, '\n');
         end
         fclose(newfile);
         
         meanfile = fopen('Mean_summary.dat', 'wt');
         for jdx = 1:numel((vis2_loc));
             fprintf(meanfile, '%f', (Means(jdx) * (180/pi)));
             fprintf(meanfile, '\n');
         end
         fclose(meanfile);
         
         varfile = fopen('Variance_summary.dat', 'wt');
         for jdx = 1:numel((vis2_loc));
             fprintf(varfile, '%f', Vars(jdx));
             fprintf(varfile, '\n');
         end
         fclose(varfile);
         
         SDfile = fopen('SD_summary.dat', 'wt');
         for jdx = 1:numel((vis2_loc));
             fprintf(SDfile, '%f', (SDs(jdx) * (180/pi)));
             fprintf(SDfile, '\n');
         end
         fclose(SDfile);
                   
         
         
pooled_RVL = (sum(ResVecLens(:)))/numel(ResVecLens);

pooled_VAR = 1-pooled_RVL;

pooled_SD = sqrt(2 * pooled_VAR);


pooledfile = fopen('PooledSD.dat', 'wt');
fprintf(pooledfile, '%f', (pooled_SD *(180/pi)));
fclose(pooledfile);

pSD(adx) = pooled_SD;

Pvectors_shorters = zeros(1,numel(vis2_loc));

for counter = 2:numel(Pvectors)
   Pvectors_shorters(counter-1) = Pvectors(counter); 
end

conflicts = [10:10:180];
undershoots = conflicts(:) - Pvectors_shorters(:);

correlations(adx) = corr2(undershoots',Means);


         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);

end


if summaries
    

figure();
semilogx(param_two, correlations, 'Linewidth', 2.0);
%plot(param_two, correlations, 'Linewidth', 2.0);
xlabel('Rotation Velocity (V)', 'Fontsize', 32);
ylabel('Correlation', 'Fontsize', 32);
set(gca, 'Fontsize', 24);
xlim([18,360]);
set(gca, 'XTick', [18,45,90,180,360]);


ylim([-0.01,1.001]);
set(gca,'YTick',[0,0.25,0.5,0.75,1]);
set(gcf,'Position', get(0,'Screensize'));   %Maximise figure to look good when saved.
set(gcf, 'PaperPositionMode', 'auto');      %Overwite tendency of 'saveas' command to resize figure back again.
box off;
saveas(gcf, 'deg_Corr', 'epsc');
close(gcf);


figure();
semilogx(param_two,(pSD * (180/pi)), 'Linewidth', 2.0);
%plot(param_two,(pSD * (180/pi)), 'Linewidth', 2.0);
xlabel('Rotation Velocity (V)', 'Fontsize', 32);
ylabel('Pooled SD', 'Fontsize', 32);
set(gca, 'Fontsize', 24);
xlim([18,360]);
set(gca, 'XTick', [18,45,90,180,360]);


ylim([-0.01,8.00]);
set(gca,'YTick',[0,2,4,6,8]);
set(gcf,'Position', get(0,'Screensize'));   %Maximise figure to look good when saved.
set(gcf, 'PaperPositionMode', 'auto');      %Overwite tendency of 'saveas' command to resize figure back again.
box off;
saveas(gcf, 'deg_SD', 'epsc');
close(gcf);

end



for adx = 1: numel(param_two)
    for jdx = 1: (numel(vis2_loc) +1);
   
        discrepancy(adx,jdx) = abs(discunder(adx,jdx) - discmean(adx,jdx));
        
    end
end

%Now work out average discrepancies for each adx value.

mean_disc = zeros(1, numel(param_two));

for adx = 1:numel (param_two)
   mean_disc(adx) = mean2(discrepancy(adx,:)); 
end

figure();
semilogx(param_two, mean_disc,'Linewidth', 2.0);
%plot(param_two, mean_disc,'Linewidth', 2.0);
xlabel('Rotation Velocity (V)', 'Fontsize', 32);
ylabel('Mean Discrepancy (deg)', 'Fontsize', 32);
set(gca, 'Fontsize', 24);
xlim([18,360]);
set(gca, 'XTick', [18,45,90,180,360]);


ylim([-0.02,14]);
set(gca,'YTick',[0:2:14]);
set(gcf,'Position', get(0,'Screensize'));   %Maximise figure to look good when saved.
set(gcf, 'PaperPositionMode', 'auto');      %Overwite tendency of 'saveas' command to resize figure back again.
box off;
saveas(gcf, 'deg_Disc', 'epsc');
close(gcf);



end
         