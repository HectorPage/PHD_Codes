function summarise_nopi(param_two,param_one, vis2_loc, timesteps, cells)
%N.B: timsteps = simulation time (in seconds) x 100
for adx = 1:numel(param_two)
    dirstringW = [num2str(param_two(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/_no_PI/parameter_search/MultiThread'];
    tier_1_path = [parentpath,'/',dirstringW];
 
    cd(tier_1_path); 
    
    for idx = 1:numel(param_one)
        dirstringS = [num2str(param_one(idx))];
        tier_2_path = [tier_1_path,'/',dirstringS];
    
        cd(tier_2_path); 
        
        %Pvectors = zeros(1,numel(vis2_loc)+1);
        %Pvectors(1) = 0;
        StandDevs = zeros(1,numel(vis2_loc));
        
        
         for jdx = 1:numel(vis2_loc)
    
            dirstringL = [num2str(vis2_loc(jdx))];
            tier_3_path = [tier_2_path, '/', dirstringL];
            cd(tier_3_path);
            
            fid = fopen('HDRates.bdat', 'rb');
            rates = fread(fid, [timesteps, cells], 'float32')';
            fclose(fid);
            
            prefdirs = linspace(1, 360, cells);
            
            numerator = 0;
            denominator = 0;
            
            Ri = zeros(2,cells);
            VecAvg = zeros(2,1);
            
            for cell = 1:cells;
                numerator = (rates(cell, timesteps)*sind(prefdirs(cell)));
                denominator = (rates(cell, timesteps)*cosd(prefdirs(cell)));
            
                Ri(1,cell) = denominator;
                Ri(2,cell) = numerator;
            end
            
            for cell = 1:cells;
                VecAvg(1,1) = VecAvg(1,1) + Ri(1,cell);
                VecAvg(2,1) = VecAvg(2,1) + Ri(2,cell);
            end
            
            VecAvg = VecAvg/500;
            
            ResVecLength = sqrt((VecAvg(1,1)*VecAvg(1,1))+(VecAvg(2,1)*VecAvg(2,1)));
                            
        
            StandDevs(jdx) = sqrt(2*(1-ResVecLength));
          
            
            cd(tier_2_path);
            
         end
         
         %x = [0:10:180];
         
         %figure();
         %plot(x,Pvectors, 'LineWidth', 2);
         %set(gca, 'FontSize',24);
         %hold on;
         %plot([0,180],[0,180],'k', 'LineWidth', 2);
         %xlabel('Conflict (Deg)', 'FontSize', 32);
         %ylabel('Packet Rotation (Deg)', 'FontSize', 32);
         %ylim([0,180]);
         %xlim([0,180]);
         %set(gca, 'FontSize',24);
         %set(gca, 'YTick',[0:60:180]);
         %set(gca, 'XTick',[0:60:180]);
         %saveas(gcf,'summary_graph', 'pdf');
         %close(gcf);
         
         
         sdfile = fopen('sd_summary.dat', 'wt');
         for jdx = 1:numel(vis2_loc)
             fprintf(sdfile, '%f', StandDevs(jdx));
             fprintf(sdfile, '\n');
         end
         fclose(sdfile);

         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);

end

end
         