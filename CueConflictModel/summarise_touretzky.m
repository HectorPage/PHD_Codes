function summarise_touretzky(light_width, RC_width, vis2_loc)

for adx = 1:numel(light_width)
    dirstringW = [num2str(light_width(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/parameter_search/sandpit/separate_rc_visring_sigmas'];
    tier_1_path = [parentpath,'/',dirstringW];
 
    cd(tier_1_path); 
    
    for idx = 1:numel(RC_width)
        dirstringS = [num2str(RC_width(idx))];
        tier_2_path = [tier_1_path,'/',dirstringS];
    
        cd(tier_2_path); 
        
        Pvectors = zeros(1,numel(vis2_loc)+1);
        Pvectors(1) = 0;
        
         for jdx = 1:numel(vis2_loc)
    
            dirstringL = [num2str(vis2_loc(jdx))];
            tier_3_path = [tier_2_path, '/', dirstringL];
            cd(tier_3_path);
            
            fid = fopen('Pvector.dat', 'r');
            str = fgets(fid);
            Pvectors(jdx+1) = sscanf(str, 'End population vector is: %f');
            Pvectors(jdx+1) = Pvectors(jdx+1) - 180;
            fclose(fid);
            
            cd(tier_2_path);
            
         end
         
         x = [0:20:180];
         
         figure();
         plot(x,Pvectors);
         hold on;
         plot([0,180],[0,180],'k');
         xlabel('Conflict (Deg)');
         ylabel('Packet Rotation (Deg)');
         ylim([0,180]);
         title(['Light Width: ',int2str(light_width(adx)),' RC Width: ', num2str(RC_width(idx))])
         saveas(gcf,'summary_graph', 'pdf');
         close(gcf);

         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);

end

end
         
         