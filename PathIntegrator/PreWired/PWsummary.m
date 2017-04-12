function PWsummary(delay, tau, symswitch)

%This script is plotting for a given target speed and given delay, what packet speed is for each strength.




packet_speed = zeros(1,numel(tau));
offset = zeros(1, numel(tau));

for adx = 1:numel(delay)
    dirstringW = [num2str(delay(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/toy_model/MultiThread/_tanh/_RC_effect/pre-wired/_symmetrical'];
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
            
            %packet_speed(jdx) = (packet_speed(jdx)/tau(jdx)) *100;
            
              fid = fopen('offset.dat', 'r');
              str = fgets(fid);
              offset(jdx) = sscanf(str, 'observed offset = %f')';
              fclose(fid);
            
            
             
            
            
            cd(tier_2_path);
            
         end
         
         
         


%          figure();
%          plot(tau, packet_speed,'Linewidth', 2);
%          xlabel('\tau^{HD}', 'Fontsize', 24);
%          ylabel('Packet Speed (^{\circ}/s)', 'Fontsize', 24);
%          xlim([0,0.1]);
%          set(gca, 'Xtick', [0:0.02:0.1]);
%          ylim([0,180]);
%          title(['Expected speed: 180^{\circ}/s'], 'Fontsize', 32);
%          set(gca, 'Fontsize', 24);
%          saveas(gcf,'PW_tau_speed', 'epsc');
%          close(gcf);
         
%            figure();
%            plot(tau, offset,'Linewidth', 2);
%            xlabel('\tau^{HD}', 'Fontsize', 24);
%            ylabel('Offset (^{\circ})', 'Fontsize', 24);
%            xlim([0,0.1]);
%            set(gca, 'Xtick', [60:60:360]);
%            ylim([0,2]);
%            set(gca, 'Fontsize', 24);
%            saveas(gcf,'PW_tau_offsets', 'epsc');
%            close(gcf);
         
         figure();
         plot(tau, ((offset/1.8)*100),'k','Linewidth', 2);
         hold on
         plot(tau, ((packet_speed/packet_speed(1))*100),'--sk','Linewidth', 2);
         xlabel('\lambda^{NO}', 'Fontsize', 24);
         ylabel('Percentage', 'Fontsize', 24);
         xlim([0,1]);
         ylim([0,100]);
         set(gca, 'Xtick', [0:0.25:1]);
         %xtick([0:0.25:1]);
         %ytick([0:25:100]);
         set(gca, 'Ytick', [0:25:100]);
         set(gca, 'Fontsize', 24);
         title(['Pre-wired'], 'Fontsize', 32);
         saveas(gcf,'paper_percentage_summary', 'epsc');
         close(gcf);

         cd(tier_1_path);
            
      
    end
    
    cd(parentpath);
    

    
    

end

end
         
