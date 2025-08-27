% script to plot


filePathOrig="./2023_04_07_ParametricStudy_Steel";

numFiles=14;
figure

p = zeros(51,1);

p100=[];
p150=[];
p200=[];
p250=[];
p300=[];

for trialPlot= 0:numFiles
    
    mat = readmatrix(filePathOrig+"/trial_"+string(round(trialPlot,0))+"/stress_strain.csv");

    p = p+mat(:,2)/3/10^6;
    
    if trialPlot>-1 && trialPlot<3
        %plot(mat(:,1),mat(:,2)/10^6,'r','LineWidth',2,'DisplayName','P=100W');
        plot(mat(:,1),mat(:,2)/10^6,'r','LineWidth',2,'DisplayName',string(trialPlot));
        hold on
        
        p100=[p100,mat(:,2)];
    end
    
    if trialPlot>2 && trialPlot<6
        if trialPlot==3
            continue
        end
        %plot(mat(:,1),mat(:,2)/10^6,'k','LineWidth',2,'DisplayName','P=150W');
        plot(mat(:,1),mat(:,2)/10^6,'k','LineWidth',2,'DisplayName',string(trialPlot));
        hold on
        
        p150=[p150,mat(:,2)];
    end
    
    if trialPlot>5 && trialPlot<9
        %plot(mat(:,1),mat(:,2)/10^6,'b','LineWidth',2,'DisplayName','P=200W');
        plot(mat(:,1),mat(:,2)/10^6,'b','LineWidth',2,'DisplayName',string(trialPlot));
        hold on
        
        p200=[p200,mat(:,2)];
    end
    
    if trialPlot>8 && trialPlot<12
        if trialPlot==11
            continue
        end
        %plot(mat(:,1),mat(:,2)/10^6,'g','LineWidth',2,'DisplayName','P=250W');
        plot(mat(:,1),mat(:,2)/10^6,'g','LineWidth',2,'DisplayName',string(trialPlot));
        hold on
        
        p250=[p250,mat(:,2)];
    end    
    
    if trialPlot>11 && trialPlot<15
        %plot(mat(:,1),mat(:,2)/10^6,'c','LineWidth',2,'DisplayName','P=250W');
        plot(mat(:,1),mat(:,2)/10^6,'c','LineWidth',2,'DisplayName',string(trialPlot));
        hold on
        
        p300=[p300,mat(:,2)];
    end   

%     if rem(trialPlot+1,3)==0
%         plot(mat(:,1),p,'LineWidth',2)
%         p = zeros(51,1);
%         hold on
% 
%     end
    
end

legend