function Post_processing(prop_save)
    close all;
    out = false; n = 1;
    % Unloads PropSim results
    t_vec = prop_save(2,2:prop_save(1,2)+1);
    P_vec = prop_save(3,2:prop_save(1,2)+1)*10^(-5);
    T_vec = prop_save(4,2:prop_save(1,2)+1);
    CF_vec = prop_save(5,2:prop_save(1,2)+1);
    It = prop_save(6,2);
    It = round(It,1);
    
    sw_results = input('Do you want to only see the results or compare them? (1 - See / 2 - Compare two plots / 3 - Compare multiple plots): ');
    switch sw_results
    case 1  % Just a direct plot of PropSim results
    figure
    
    ax1 = subplot(2,3,1); % Left    
    plot(t_vec,P_vec);
    title('Chamber pressure');
    xlabel(['Time (s)']);
    ylabel(['Chamber pressure (bar)']);
    
    ax2 = subplot(2,3,2); % Middle
    plot(t_vec,T_vec);
    title('Thrust');
    xlabel(['Time (s)']);
    ylabel(['Thrust (N)']); 
    
    ax3 = subplot(2,3,3); % Right  
    plot(t_vec,CF_vec);
    title('Thrust coefficient');
    xlabel(['Time (s)']);
    ylabel(['Thrust coefficient']);
    
    sgtitle(['PropSim simulation results (It = ',num2str(It),' Ns)'])
    set(gcf,'WindowState','maximized')
    
    case 3
        
    ax1 = subplot(2,2,1); % Left upper
    plot(t_vec,T_vec);
    title('Thrust');
    xlabel(['Time (s)']);
    ylabel(['Thrust (N)']); 
    
    ax2 = subplot(2,2,2); % Right upper  
    plot(t_vec,CF_vec);
    title('Thrust coefficient');
    xlabel(['Time (s)']);
    ylabel(['Thrust coefficient']);
    
    ax3 = subplot(2,2,3); % Left lower   
    plot(t_vec,P_vec,'DisplayName',['PropSim: It = ',num2str(It),' Ns']);
    title('Chamber pressure');
    xlabel(['Time (s)']);
    ylabel(['Chamber pressure (bar)']);
    
    set(gcf,'WindowState','maximized')
    set(gcf,'color','w'); 
    
    while ~out % Acumulate plots
    fprintf('Choose a csv curve with which to compare the results.\n')
    [file,path] = uigetfile('.\Saved curves\*.csv');
    fullpath = cat(2,path,file);
    prop_comp = readtable(fullpath);
    %test = prop_comp{:,2:prop_comp{1,2}};
    t_vec_comp = prop_comp{2,2:prop_comp{1,2}+1};
    if prop_comp{3,2} > 1000 % Pressure must be in Pa
        P_vec_comp = prop_comp{3,2:prop_comp{1,2}+1}*10^(-5);
    else % Pressure must be in bar
        P_vec_comp = prop_comp{3,2:prop_comp{1,2}+1};
    end
    T_vec_comp = prop_comp{4,2:prop_comp{1,2}+1};
    CF_vec_comp = prop_comp{5,2:prop_comp{1,2}+1};
    It_comp = prop_comp{6,2};
    It_comp = round(It_comp,1);     
    
    hold(ax1,'on')
    plot(ax1,t_vec_comp,T_vec_comp);

    hold(ax2,'on')
    plot(ax2,t_vec_comp,CF_vec_comp);
    
    hold(ax3,'on')
    plot(ax3,t_vec_comp,P_vec_comp,'DisplayName',[file(1:end-4),': It = ',num2str(It_comp),' Ns']);   
    
    legend('Position',[0.6 0.2 0.2 0.2]) % Legend on the right lower
    
    new_comp = input('Do you want to compare against a different graph? (Y/N): ','s');
    if new_comp == 'y' || new_comp == 'Y'
        out = false;       
    else
        out = true;
    end 
    end
    
    otherwise

    while out == false
    fprintf('Choose a csv curve with which to compare the results.\n')
    [file,path] = uigetfile('.\Saved curves\*.csv');
    fullpath = cat(2,path,file);
    prop_comp = readtable(fullpath);
    t_vec_comp = prop_comp{2,2:prop_comp{1,2}+1};
    if prop_comp{3,2} > 1000 % Pressure must be in Pa
        P_vec_comp = prop_comp{3,2:prop_comp{1,2}+1}*10^(-5);
    else % Pressure must be in bar
        P_vec_comp = prop_comp{3,2:prop_comp{1,2}+1};
    end
    
    T_vec_comp = prop_comp{4,2:prop_comp{1,2}+1};
    CF_vec_comp = prop_comp{5,2:prop_comp{1,2}+1};
    It_comp = prop_comp{6,2};
    It_comp = round(It_comp,1);    
    
    rel_error = abs(It-It_comp)/It_comp*100;
    
    figure(n)   
    ax1 = subplot(3,2,1); % Left upper
    plot(t_vec,P_vec,'');
    hold(ax1,'on')
    plot(ax1,t_vec_comp,P_vec_comp);
    legend('PropSim',file(1:end-4))
    xlabel(['Time (s)']);
    ylabel(['Chamber pressure (bar)']);

    ax3 = subplot(3,2,3); % Left middle
    plot(t_vec,T_vec);
    hold(ax3,'on')
    plot(ax3,t_vec_comp,T_vec_comp);
    legend('PropSim',file(1:end-4))
    xlabel(['Time (s)']);
    ylabel(['Thrust (N)']);
    
    
    ax5 = subplot(3,2,5); % Left lower
    plot(t_vec,CF_vec);
    hold(ax5,'on')
    plot(ax5,t_vec_comp,CF_vec_comp);
    legend('PropSim',file(1:end-4))
    xlabel(['Time (s)']);
    ylabel(['Thrust coefficient']);
    
    P_vec_comp = P_vec_comp(P_vec_comp ~= 0);
    T_vec_comp = T_vec_comp(1:length(P_vec_comp));
    CF_vec_comp = CF_vec_comp(1:length(P_vec_comp));
    t_vec_comp = t_vec_comp(1:length(P_vec_comp));
    
    if t_vec(end) > t_vec_comp(end)
        % Simulation curve longer than reference curve
        t_vec_corr = t_vec_comp;
        for i = 1:length(t_vec_corr)
            % You correct the curve given by PropSim to follow reference
            P_vec_corr(i) = interp1(t_vec,P_vec,t_vec_corr(i));
            T_vec_corr(i) = interp1(t_vec,T_vec,t_vec_corr(i)); 
            CF_vec_corr(i) = interp1(t_vec,CF_vec,t_vec_corr(i)); 
            
            err_P(i) = abs(P_vec_corr(i) - P_vec_comp(i))/P_vec_comp(i)*100;
            err_T(i) = abs(T_vec_corr(i) - T_vec_comp(i))/T_vec_comp(i)*100;
            err_CF(i) = abs(CF_vec_corr(i) - CF_vec_comp(i))/CF_vec_comp(i)*100;
            
            if err_P(i) > 100 || isnan(err_P(i))
                err_P(i) = 100;
            elseif err_T(i) > 100 || isnan(err_T(i))
                err_T(i) = 100;
            elseif err_CF(i) > 100 || isnan(err_CF(i))
                err_CF(i) = 100;
            else
            end
        end   
    else % Reference curve longer than simulation curve
        t_vec_corr = t_vec;
        for i = 1:length(t_vec_corr)
            % You correct the reference to follow PropSim
            P_vec_corr(i) = interp1(t_vec_comp,P_vec_comp,t_vec_corr(i)); % Resizing sim curve as reference
            T_vec_corr(i) = interp1(t_vec_comp,T_vec_comp,t_vec_corr(i));
            CF_vec_corr(i) = interp1(t_vec_comp,CF_vec_comp,t_vec_corr(i));
            
            err_P(i) = abs(P_vec_corr(i) - P_vec(i))/P_vec(i)*100;
            err_T(i) = abs(T_vec_corr(i) - T_vec(i))/T_vec(i)*100;
            err_CF(i) = abs(CF_vec_corr(i) - CF_vec(i))/CF_vec(i)*100;
            
            if err_P(i) > 100 || isnan(err_P(i))
                err_P(i) = 100;
            elseif err_T(i) > 100 || isnan(err_T(i))
                err_T(i) = 100;
            elseif err_CF(i) > 100 || isnan(err_CF(i))
                err_CF(i) = 100;
            else
            end
        end
    end

    ax2 = subplot(3,2,2); % Right upper
    plot(t_vec_corr,err_P);
    xlabel(['Time (s)']);
    ylabel(['Relative error pressure (%)']);
    ylim([0 100]);
    xlim([0 t_vec_corr(end)+0.5]);
    
    ax4 = subplot(3,2,4); % Right middle
    plot(t_vec_corr,err_T);
    xlabel(['Time (s)']);
    ylabel(['Relative error thrust (%)']);
    ylim([0 100]);
    xlim([0 t_vec_corr(end)+0.5]);
    
    ax6 = subplot(3,2,6); % Right lower
    plot(t_vec_corr,err_CF);
    xlabel(['Time (s)']);
    ylabel(['Relative error thrust coefficient (%)']);
    ylim([0 100]);
    xlim([0 t_vec_corr(end)+0.5]);
    
    sgtitle(['PropSim results (It = ',num2str(It),') compared to ',file(1:end-4),' curve (It = ',num2str(It_comp),') (relative error It = ',num2str(rel_error),'%)'])
    set(gcf,'WindowState','maximized')
    set(gcf,'color','w')
    
    new_comp = input('Do you want to compare against a different graph? (Y/N): ','s');
    if new_comp == 'y' || new_comp == 'Y'
        out = false;
        new_fig = input('Do you want to keep the previous plot? (Y/N): ','s');
        if new_fig == 'y' || new_fig == 'Y'
            n = n+1;
        else
            close all;
        end
        
    else
        out = true;
    end
    end
       

    end
end