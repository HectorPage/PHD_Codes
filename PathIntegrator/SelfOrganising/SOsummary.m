function SOsummary(delay, tau, symswitch)

%This script is plotting for a given target speed and given delay, what packet speed is for each strength.




packet_speed = zeros(1,numel(tau));
offset = zeros(1, numel(tau));

for adx = 1:numel(delay)
    dirstringW = [num2str(delay(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/toy_model/MultiThread/_tanh/_RC_effect/nonrandom'];
    tier_1_path = [parentpath,'/',dirstringW];
 
    cd(tier_1_path); 
    
    for idx = 1:numel(symswitch)
        dirstringS = [num2str(symswitch(idx))];
        tier_2_path = [tier_1_path,'/',dirstringS];
    
        cd(tier_2_path); 
        
        
         for jdx = 1:numel(tau)
             
            
    
            dirstringL = [num2str(tau(jdx))];
            tier_3_path = [tier_2_path, '/sym', dirstringL];
            cd(tier_3_path);
            
            %FOLLOWING PIECE OF CODE HERE TO EXTRACT PACKET SPEED AND OFFSET.
            
            
            fid = fopen('speed.dat', 'r');
            str = fgets(fid);
            packet_speed(jdx) = sscanf(str, 'speed: %f');
            fclose(fid);
            
            %packet_speed(jdx) = (packet_speed(jdx)/delay(jdx)) *100;
            %THIS LINE ONLY USED WHEN VARYING TARGET SPEED!
            
              fid = fopen('offset_mean.dat', 'r');
              str = fgets(fid);
              offset(jdx) = sscanf(str, '%f')';
              fclose(fid);
            
            
   
            cd(tier_2_path);
            
            %ONLY FOR THE SPEED PLOTS
            %packet_speed(jdx) = packet_speed(jdx)/tau(jdx) * 100;
            
         end
         
         
         

% 
%          figure();
%          plot(tau, packet_speed ,'Linewidth', 2);
%          xlabel('\lambda^{NO}', 'Fontsize', 24);
%          ylabel('Packet Speed (^{\circ}/s) ', 'Fontsize', 24);
%          xlim([0,1]);
%          set(gca, 'Xtick', [0:0.25:1]);
%          ylim([0,180]);
%          title(['Expected speed: 180^{\circ}/s'], 'Fontsize', 32);
%          set(gca, 'Fontsize', 24);
%          saveas(gcf,'SO_sym_speed', 'epsc');
%          close(gcf);
%          
%            figure();
%             plot(tau, offset,'Linewidth', 2);
%           xlabel('\lambda^{NO}', 'Fontsize', 24);
%             ylabel('Offset (^{\circ})', 'Fontsize', 24);
%           xlim([0,1]);
%        set(gca, 'Xtick', [0:0.25:1]);
%             ylim([0,1.8]);
%             set(gca, 'Fontsize', 24);
%             title(['Expected offset: 1.8^{\circ}'], 'Fontsize', 32);
%             saveas(gcf,'SO_sym_offset', 'epsc');
%             close(gcf);
%          
           figure();
           plot(tau,((offset/offset(1))*100),'k','Linewidth', 2);
           hold on
           plot(tau, ((packet_speed/packet_speed(1))*100),'--sk','Linewidth', 2);
           xlabel('\lambda^{NO}', 'Fontsize', 24);
           ylabel('Percentage', 'Fontsize', 24);
            xlim([0,1]);
            set(gca, 'Xtick', [0:0.25:1]);
           ylim([0,100]);
           set(gca, 'Fontsize', 24);
           title(['Expected offset: 1.8^{\circ}'], 'Fontsize', 32);
           saveas(gcf,'SO_sym_percentage', 'epsc');
           close(gcf);

         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);
    

    
    

end

end
         
