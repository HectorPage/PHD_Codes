function NORC_PWsummary(delay, tau, symswitch)

%This script is plotting for a given target speed and given delay, what packet speed is for each strength.




packet_speed = zeros(1,numel(tau));
phase_offset = zeros(1,numel(tau));
offset = zeros(1, numel(tau));

for adx = 1:numel(delay)
    dirstringW = [num2str(delay(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/no_RC/single_layer/MultiThread/BothDelays/parameter_search'];
    tier_1_path = [parentpath,'/',dirstringW];
 
    cd(tier_1_path); 
    
    for idx = 1:numel(symswitch)
        dirstringS = [num2str(symswitch(idx))];
        tier_2_path = [tier_1_path,'/',dirstringS];
    
        cd(tier_2_path); 
        
        
         for jdx = 1:numel(tau)
             
            
    
            dirstringL = [num2str(tau(jdx))];
            tier_3_path = [tier_2_path, '/delay', dirstringL];
            cd(tier_3_path);
            
            %FOLLOWING PIECE OF CODE HERE TO EXTRACT PACKET SPEED AND OFFSET.
            
            
            fid = fopen('speed.dat', 'r');
            str = fgets(fid);
            packet_speed(jdx) = sscanf(str, 'speed: %f');
            fclose(fid);
            
%             fid = fopen('_PhaseOffset.dat', 'r');
%             str = fgets(fid);
%             phase_offset(jdx) = sscanf(str, '%f');
%             fclose(fid);
            
            
            %packet_speed(jdx) = (packet_speed(jdx)/tau(jdx)) *100;
            
%               fid = fopen('offset.dat', 'r');
%               str = fgets(fid);
%               offset(jdx) = sscanf(str, 'observed offset = %f')';
%               fclose(fid);
%             
            
             
            
            
            cd(tier_2_path);
            
         end
         
         
         


%          figure();
%          plot(tau, packet_speed,'Linewidth', 2);
%          xlabel('Target Velocity', 'Fontsize', 24);
%          ylabel('Precentage', 'Fontsize', 24);
%          xlim([45,360]);
%          set(gca, 'Xtick', [0:90:360]);
%          ylim([0,110]);
%          set(gca, 'Ytick', [0:25:100]);
%          title('Observed Velocity as % of Target', 'Fontsize', 32);
%          set(gca, 'Fontsize', 24);
%          saveas(gcf,'PW_velocity_speed', 'epsc');
%          close(gcf);
         
           figure();
           plot(tau, packet_speed,'Linewidth', 2);
           xlabel('\Delta t', 'Fontsize', 24);
           ylabel('Speed (deg/s)', 'Fontsize', 24);
           xlim([0,0.05]);
           ylim([0,180]);
           set(gca, 'Ytick', [0:20:180]);
           set(gca, 'Fontsize', 24);
           saveas(gcf,'PW_delay_speeds', 'epsc');
           close(gcf);
%          
%          figure();
% %        plot(tau, ((offset/1.8)*100),'k','Linewidth', 2);
% %        hold on
%          plot(tau, packet_speed, 'Linewidth', 2);
%          xlabel('Conduction Delay (s)', 'Fontsize', 24);
%          ylabel('Speed (deg/s)', 'Fontsize', 24);
%          set(gca, 'Fontsize', 24);
%          title(['Path Integration Speed'], 'Fontsize', 32);
%          saveas(gcf,'delay_speed_summary', 'epsc');
%          close(gcf);
%          
%          
%          figure();
%          plot(tau, phase_offset, 'Linewidth',2);
%          %ylim([0, 0.21]);
%          %set(gca, 'Ytick',[0.025, 0.05, 0.075, 0.1]);
%          %xlim([0, 0.11]);
%          %set(gca, 'Xtick',[0.025, 0.05, 0.075, 0.1]);
%          xlabel('Conduction Delay (s)', 'Fontsize', 24);
%          ylabel('Phase Offset (s)', 'Fontsize', 24);
%          set(gca, 'Fontsize', 24);
%          title(['Phase Offset'], 'Fontsize', 32);
%          saveas(gcf,'delay_phase_summary', 'epsc');
%          close(gcf);

         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);
    

    
    

end

end
         
