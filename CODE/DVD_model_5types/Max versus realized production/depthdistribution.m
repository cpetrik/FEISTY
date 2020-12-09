function depthdistribution(param)
    % calculate depth distribution per meter
    xrange = linspace(0,param.bottom,(param.bottom+1));
    
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
    if ~make_it_tight,  clear subplot;  end

    figure(3)
    % zooplankton
    subplot(9,6,1)
    barh(-xrange,param.depthDay(1:end,1),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.02])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,2)
    barh(-xrange,param.depthDay(1:end,2),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,2),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.02])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    % benthos
    subplot(9,6,7)
    barh(-xrange,param.depthDay(1:end,3),'FaceColor',[0.494 0.184 0.556]) 
    xticks([0 0.1])
    ylim([-max(xrange) 0])
    yticks([-max(xrange) 0])  
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,8)
    barh(-xrange,param.depthDay(1:end,4),'FaceColor',[0.494 0.184 0.556]) 
    xticks([0 0.1])
    ylim([-max(xrange) 0])
    yticks([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    % small pelagics
    for i = 1:4
        subplot(9,6,12+i)
        barh(-xrange,param.depthDay(1:end,4+i),'FaceColor',[0.494 0.184 0.556]) 
        xlim([0 0.01])
        ylim([-max(xrange) 0])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end
      
    % meso pelagics
    subplot(9,6,19)
    barh(-xrange,param.depthDay(1:end,9),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
        
    subplot(9,6,20)
    barh(-xrange,param.depthDay(1:end,10),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[]) 
    
    subplot(9,6,21)
    barh(-xrange,param.depthDay(1:end,11),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,11),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[]) 
   
    subplot(9,6,22)
    barh(-xrange,param.depthDay(1:end,12),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,12),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    % large pelagics
    for i = 1:6
        subplot(9,6,24+i)
        barh(-xrange,param.depthDay(1:end,12+i),'FaceColor',[0.494 0.184 0.556]) 
        xlim([0 0.01])
        ylim([-max(xrange) 0])
        set(gca,'YTick',[])
        set(gca,'XTick',[]) 
    end
      
    % large bathypelagics
    subplot(9,6,31)
    barh(-xrange,param.depthDay(1:end,19),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[]) 
    
    subplot(9,6,32)
    barh(-xrange,param.depthDay(1:end,20),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,33)
    barh(-xrange,param.depthDay(1:end,21),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,21),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,34)
    barh(-xrange,param.depthDay(1:end,22),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,22),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,35)
    barh(-xrange,param.depthDay(1:end,23),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,36)
    barh(-xrange,param.depthDay(1:end,24),'FaceColor',[0.494 0.184 0.556]) 
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    % large demersal 
     for i = 1:4
        subplot(9,6,36+i)
        barh(-xrange,param.depthDay(1:end,24+i),'FaceColor',[0.494 0.184 0.556]) 
        xlim([0 0.01])
        ylim([-max(xrange) 0])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end
    
    subplot(9,6,41)
     barh(-xrange,param.depthDay(1:end,29),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,29),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(9,6,42)
    barh(-xrange,param.depthDay(1:end,30),'FaceColor',[0 0.447 0.741]) 
    hold on
    barh(-xrange,param.depthNight(1:end,30),'FaceColor',[0.85 0.325 0.098])
    hold off
    xlim([0 0.01])
    ylim([-max(xrange) 0])
    set(gca,'YTick',[])
    set(gca,'XTick',[])