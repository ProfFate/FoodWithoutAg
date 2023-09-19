function [Check] = FoodWithoutAg()

% Code to generate plots for Davis, SJ, K Alexander, J Moreno-Cruz, C Hong,
% M Shaner, K Caldeira, and I McKay. (2023) Food Without Agriculture.
% Nature Sustainability, v. x, p. xxx-xxx

clear

%% Initialize Variables
PrintFigs = 1;

n = 1e2;

% Emissions intensity of energy
minEmsIntensity_energy = 0; % gCO2/MJ
maxEmsIntensity_energy = 282; % gCO2/MJ

%Agriculture

    % Land-use emissions
    minLUEms = 0.25; % oil crops ROW from Hong et al. 2021, gCO2e/kcal.
    maxLUEms = 2.5; % oil crops sub-Saharan Africa from Hong et al. 2021, gCO2e/kcal.

    % Energy intensity for ag
    minEngyIntensity_ag = 0.0006; %MJ/kcal  i.e. low end for soy from Crippa et al 2021
    maxEngyIntensity_ag = 0.0102; %MJ/kcal  i.e. high end for palm from Crippa et al 2021
    
% Synthetic
    % Feedstock emissions (vertical in b)
    minEmsIntensity_feed = 0; %i.e. CO2. kgCO2/kg feed (sum of emissions values in Table S7)
    midEmsIntensity_feed = 3.16; %i.e. gas. kgCO2/kg feed (sum of emissions values in Table S7)
    maxEmsIntensity_feed = 3.51; %i.e. coal. kgCO2/kg feed (sum of emissions values in Table S7)

    minFeedstockDemand = 4.07; %i.e. CO2. kg feed/kg product (inverse of eta_feedstock values in Table S7)
    midFeedstockDemand = 1.51; %i.e. gas. kg feed/kg product (inverse of eta_feedstock values in Table S7)
    maxFeedstockDemand = 2.91; %i.e. coal. kg feed/kg product (inverse of eta_feedstock values in Table S7)
    
    % Energy emissions (horizontal in b)
    minEngyIntensity_chem = 4.9; %i.e. coal. kWh/kg product (converted from production energy values in Table S7)
    midEngyIntensity_chem = 4.8; %i.e. gas. kWh/kg product (converted from production energy values in Table S7)
    maxEngyIntensity_chem = 39.3; %i.e. DAC. kWh/kg product (converted from production energy values in Table S7)

%% Calculate emissions intensities as functions of energy emissions and input emissions

% Agriculture
range_LUEms = maxLUEms - minLUEms;
minEngyEms_ag = round(minEmsIntensity_energy * minEngyIntensity_ag,1);
%maxEngyEms_ag = round(maxEmsIntensity_energy * maxEngyIntensity_ag,1);
maxEngyEms_ag = 1.5; 

range_EngyEms_ag = maxEngyEms_ag - minEngyEms_ag;
increment_LUEms = range_LUEms/n;
increment_EngyEms_ag = range_EngyEms_ag/n;

AgEmsIntensity = [];
for i=1:n
    for j=1:n
        AgEmsIntensity(i,j) = (minLUEms + (increment_LUEms * (i-1))) + (minEngyEms_ag + (increment_EngyEms_ag * (j-1)));
    end
end

% Chemistry
FeedstockIntensity = [];
EngyIntensity = []; %kWh per kg
k2k = 0.11111; %g prod per kcal (assumed energy density of produced fat)
for j=1:n
    if j==1
       FeedstockIntensity(j) = minEmsIntensity_feed * minFeedstockDemand * k2k; % gCO2/kcal = (kgCO2/kg feed) * (kg feed/kg product) * (g prod per kcal)
       EngyIntensity(j) = minEngyIntensity_chem; % kWh/kg product
    else if j<(n/2)
            FeedstockIntensity(j) = FeedstockIntensity(j-1) + ...
                          ( ((midEmsIntensity_feed * midFeedstockDemand * k2k) - (minEmsIntensity_feed * minFeedstockDemand * k2k)) /((n/2)-1) );
            EngyIntensity(j) = minEngyIntensity_chem;
         else if j==(n/2)
                 FeedstockIntensity(j) = midEmsIntensity_feed * midFeedstockDemand * k2k;
                 EngyIntensity(j) = midEngyIntensity_chem;
              else if j>(n/2)
                      FeedstockIntensity(j) = FeedstockIntensity(j-1) + ...
                        (  ((maxEmsIntensity_feed * maxFeedstockDemand * k2k) - (midEmsIntensity_feed * midFeedstockDemand * k2k)) /(n/2) );
                      EngyIntensity(j) = EngyIntensity(j-1) + (maxEngyIntensity_chem - midEngyIntensity_chem)/(n/2);
                  end
             end
         end
    end
end
EngyIntensity = fliplr(EngyIntensity);


EngyIntensity = EngyIntensity .* 3.6 ./ 9000; %convert from kWh/kg to MJ/kcal
minEngyEms_chem = round(minEmsIntensity_energy * min(EngyIntensity),1); % gCO2/kcal = (gCO2/MJ) * (MJ/kcal)
maxEngyEms_chem = round(maxEmsIntensity_energy * max(EngyIntensity),1);

range_EmsIntensity_energy = maxEmsIntensity_energy - minEmsIntensity_energy;
increment_EmsIntensity_energy = range_EmsIntensity_energy/(n-1);


ChemEmsIntensity = [];
for i=1:n
    for j=1:n
       ChemEmsIntensity(i,j) = FeedstockIntensity(i) + ( EngyIntensity(i) * (minEmsIntensity_energy + (increment_EmsIntensity_energy * (j-1))) );
    end
end

 maxFeedstockEms = 1.2;
 
%% Set Values of Points

AgLit = []; %col 1 is ems intensity of engy, col 2 is LUems, col 3 is energy intensity of production, col4 is crop type (1 is general, 2 is palm, 3 is soy)
oilEF = 218; %gCO2/MJ, incl. life cycle

%global
AgLit(1,:) = [oilEF 0.56984 0.0005 1];
% gCO2/MJ LCA oil ems intensity
% gCO2e/kcal Europe oilcrop avg from Hong et al 2021
% MJ/kcal dervied from dividing 0.041 gCO2/kcal (low energy from Crippa et al 2021) by assumed ems intensity

AgLit(end+1,:) = [oilEF 2.28 0.0013 1];
% gCO2/MJ LCA oil ems intensity 
% gCO2e/kcal SSA from Hong et al 2021
% MJ/kcal dervied from dividing 0.1 gCO2/kcal (high energy from Crippa et al 2021) by assumed ems intensity


% Palm Oil
    %Brazil
    AgLit(end+1,:) = [oilEF 0.2425 0.0053 2];
    % gCO2/MJ LCA oil ems intensity 
    % gCO2e/kcal Brazil oil palm from Hong et al 2021
    % avg MJ/kcal from 3 farms in Angarita et al 2009

    %Colombia
    AgLit(end+1,:) = [oilEF 0.2155 0.003 2];
    % gCO2/MJ LCA oil ems intensity 
    % gCO2e/kcal Colombia palm oil from Hong et al 2021
    % avg MJ/kcal from 2 farms in Angarita et al 2009
  
    %Indonesia
    AgLit(end+1,:) = [oilEF 0.3566 0.0053 2];
    % gCO2/MJ LCA oil ems intensity 
    % gCO2e/kcal Indonesia palm oil from Hong et al 2021
    % avg MJ/kcal from Kamahara et al. 2010

    %SSA
    AgLit(end+1,:) = [oilEF 2.28 0.0038 2];
    % gCO2/MJ LCA oil ems intensity 
    % gCO2e/kcal SSA from Hong et al 2021
    % MJ/kcal dervied from Achten et al LCA of Camaroonian oil palm by assumed ems intensity

%     Malaysia
%     AgLit(end+1,:) = [oilEF 0.1440 0.003 2];
%     gCO2/MJ LCA oil ems intensity
%     gCO2e/kcal Malaysia palm oil from Hong et al 2021
%     avg MJ/kcal from *could not find*
     
% Soybean Oil
    % Europe
    AgLit(end+1,:) = [oilEF 0.56984 0.0007 3];
    % gCO2/MJ LCA oil ems intensity
    % gCO2e/kcal Europe oilcrop avg from Hong et al 2021
    % avg MJ/kcal average from Alluvione et al 2011  

    % US
    AgLit(end+1,:) = [oilEF 0.32225 0.0008 3];
    % gCO2/MJ LCA oil ems intensity
    % gCO2e/kcal US soy from Hong et al 2021
    % avg MJ/kcal average from Pradhan et al 2011  
    
    % Brazil
    AgLit(end+1,:) = [oilEF 1.574 0.0009 3];
    % gCO2/MJ LCA oil ems intensity
    % gCO2e/kcal Brazil soy from Hong et al 2021
    % avg MJ/kcal average from Zortea et al 2018
    
    % China
    AgLit(end+1,:) = [oilEF 1.04 0.00084 3];
    % gCO2/MJ LCA oil ems intensity
    % gCO2e/kcal Chinese soy from Hong et al 2021
    % avg MJ/kcal average from Hu et al 2008

SynthLit = []; %col 1 is engy ems, col 2 is feedstock ems, col 3 is product type (1 is historical fat, 2 is hypothetical fat, 3 is is protein)
               %col 4 is co2/MJ
               
    % coal feedstock, coal power
    SynthLit(1,:) = [0.77 1.13 1 278];
    
    % coal feedstock, 2020 grid mix
    SynthLit(end+1,:) = [0.77 1.13 1 111];
    
    % US: natgas feedstock, 2020 grid mix (hypothetical)
    %SynthLit(end+1,:) = [0.31 0.53 2 111];
    SynthLit(end+1,:) = [0.31 0.6 2 111];
    
    % US: CO2 feedstock, 2020 grid mix (hypothetical)
    %SynthLit(end+1,:) = [1.78 0 2 111];
    SynthLit(end+1,:) = [0 0.6 2 111];
    
    % US: natgas feedstock, decarbonized grid (hypothetical)
    SynthLit(end+1,:) = [0 0.53 2 0];
    
    % US: CO2 feedstock, decarbonized grid (hypothetical)
    SynthLit(end+1,:) = [0 0 2 0];
    
    % Methionine
    SynthLit(end+1,:) = [0 0 3 0];
    
%     SynthLit(1,1,1) = 0.041;  %methionine
%     SynthLit(1,2,1) = 1.35; % methionine, Blonk report
% 
%     SynthLit(1,1,2) = 0.041;  %engy ems, methanotroph protein, Carbon Trust
%     SynthLit(1,2,2) = 1.45; %feedstock ems, methanotroph protein, Carbon Trust

% make x,y points for SynthLit
    for i=1:size(SynthLit,2)
        SynthLit_plot(i,1) =  SynthLit(i,4) / maxEmsIntensity_energy * n;
        SynthLit_plot(i,2) = SynthLit(i,2) / maxFeedstockEms * n;
    end
    % gCO2/kcal = (gCO2/kcal) + (MJ/kcal * gCO2/MJ)
    
    
%% Figure 2 - Make contour plots of emissions/kcal

if PrintFigs==1
    
    figure

    %conventional ag oils (panel a)
    subplot(2,1,1)
%     contourf(AgEmsIntensity)
    s = pcolor(AgEmsIntensity);
    set(s, 'EdgeColor', 'none')
    s.FaceColor = 'interp';
    hold on
    contour(AgEmsIntensity, [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5], '-k');
    
    xlabel('Energy emissions - gCO2/kcal');
    ylabel('Land-use emissions - gCO2e/kcal');
    title('Agricultural oils')
    xticks([0 (n/5*1) (n/5*2) (n/5*3) (n/5*4) n])
    xticklabels({minEngyEms_ag,(maxEngyEms_ag/5*1),(maxEngyEms_ag/5*2),(maxEngyEms_ag/5*3),(maxEngyEms_ag/5*4),maxEngyEms_ag})    
    yticks([0 (n/5*1) (n/5*2) (n/5*3) (n/5*4) n])
    yticklabels({minLUEms,(maxLUEms/5*1),(maxLUEms/5*2),(maxLUEms/5*3),(maxLUEms/5*4),maxLUEms})
    
    % plot points from literature
    AgLit(:,5) = AgLit(:,1) .* AgLit(:,3); % generate field of emissions / kcal
    for i=1:size(AgLit,1)
       if AgLit(i,4)==1 %general
           hold on
           scatter((AgLit(i,5)*n/maxEngyEms_ag),(AgLit(i,2)*n/maxLUEms),'MarkerFaceColor','black','MarkerEdgeColor','none');
       else if AgLit(i,4)==2 %palm
           hold on
           scatter((AgLit(i,5)*n/maxEngyEms_ag),(AgLit(i,2)*n/maxLUEms),'MarkerFaceColor','green','MarkerEdgeColor','none');
           else if AgLit(i,4)==3 %soy
               hold on
               scatter((AgLit(i,5)*n/maxEngyEms_ag),(AgLit(i,2)*n/maxLUEms),'MarkerFaceColor','red','MarkerEdgeColor','none');
               end
           end
       end
    end %for i
    
    %Set and show colorbar
    map = flipud(cbrewer2('Spectral'));
    colormap(map)
    if max(max(AgEmsIntensity)) > max(max(ChemEmsIntensity))
        upperbound = max(max(AgEmsIntensity));
    else
        upperbound = max(max(ChemEmsIntensity));
    end
    caxis([0 round(upperbound,1)])
    colorbar %gCOe2/kcal
    
    
%synthetic oils (panel b)

    subplot(2,1,2)
    s = pcolor(ChemEmsIntensity);
    set(s, 'EdgeColor', 'none')
    s.FaceColor = 'interp';
    hold on
    contour(ChemEmsIntensity, [0.5 1 1.5 2 2.5 3 3.5 4.0], '-k');
    
    xlabel('Emissions intensity of energy - gCO2/MJ');
    ylabel('Feedstock emissions - gCO2e/kcal');
    title('Synthetic oils')
    xticks([0 (n/6*1) (n/6*2) (n/6*3) (n/6*4) (n/6*5) n])
    xticklabels({minEmsIntensity_energy,(maxEmsIntensity_energy/6*1),(maxEmsIntensity_energy/6*2),(maxEmsIntensity_energy/6*3),(maxEmsIntensity_energy/6*4),(maxEmsIntensity_energy/6*5),maxEmsIntensity_energy})
    yticks([0 (n/5*1) (n/5*2) (n/5*3) (n/5*4) n])
    yticklabels({0,(maxFeedstockEms/5*1),(maxFeedstockEms/5*2),(maxFeedstockEms/5*3),(maxFeedstockEms/5*4),maxFeedstockEms})
   
    for i=1:size(SynthLit_plot,1)
       hold on
       scatter(SynthLit_plot(i,1),SynthLit_plot(i,2),'MarkerFaceColor','black','MarkerEdgeColor','none');
    end

    
    %Set and show colorbar
    map = flipud(cbrewer2('Spectral'));
    colormap(map)
    caxis([0 upperbound])
    colorbar %gCOe2/kcal

    if PrintFigs == 1
       print(gcf,'-depsc','-painters', '../Plots/Figure2[Contours]');
    end

end %if Fig 2

%% Figure 3 - Make cumulative plots of avoided emissions and spared land 

if PrintFigs==1
    
    figure

    load CountryData.mat 
    % load data from Hong et al. 2021, including total emissions (GtCO2e), estimated avoided emissions (GtCO2e), production (kcal), and land areas (hectares)
    % by years (2010-2017 in columns), country (rows), and crops :,:,1 = palm and :,:,2 = soy

    Results = [];
    mylabels={};
    %process for table/plotting
    for cropchoice = 1:2
        Temp=[];
        tCnty = Countries(:,:,cropchoice);
        tEms = AvoidedEms(:,:,cropchoice);
        tEms = mean(tEms,2);
        tProd = Prod(:,:,cropchoice);
        tProd = mean(tProd,2);
        tArea = Area(:,:,cropchoice);
        tArea = mean(tArea,2);

        T_Results = table(tCnty,tEms,tProd,tArea);
        

        T_Results = sortrows(T_Results,'tEms','descend');
        Templabels = T_Results.tCnty;

        Temp(:,1)=T_Results.tEms;
        Temp(:,2)=T_Results.tProd;
        Temp(:,3)=T_Results.tArea;

        %cut off zeros
        i = 0;
        counter = 1;
        while i==0
            if Temp(counter,1)~=0
                i = 0;
                counter = counter + 1;
            else
                Temp = Temp(1:counter-1,:);
                Templabels = Templabels(1:counter-1,:);
                i = 1;
            end
        end

        if cropchoice==1
            CT=[0.8500, 0.3250, 0.0980]; %palm is orange
        else
            CT=[0,0.5,0]; %soy is green
        end
        barColorMap = repmat(CT,size(Temp,1),1);
        Temp(:,4:6)=barColorMap; %apply colormap

        Results = [Results; Temp];
        mylabels = {mylabels; Templabels};

    end %cropchoice

    %Panel A - avoided emissions

        % make CDF sorted by highest to lowest AvoidedEms
        Sorted = sortrows(Results,1,'descend');
        xvar = 2; %Production
        yvar = 1; %AvoidedEms
        for i = 1:size(Results,1)
            if i==1
                CDFx(i,1) = Sorted(i,xvar);
                CDFy(i,1) = Sorted(i,yvar);
            else
                CDFx(i,1) = Sorted(i,xvar) + CDFx(i-1,1);
                CDFy(i,1) = Sorted(i,yvar) + CDFy(i-1,1);
            end
        end

        subplot(2,1,1)
        edges = zeros(length(CDFx)+1,1);
        edges(2:end,1) = CDFx;
        vals = CDFy;
        center = (edges(1:end-1) + edges(2:end))/2;
        width = diff(edges);
        hold on
        for i=1:length(center)
            bar(center(i),vals(i),width(i),'FaceColor',Sorted(i,4:6),'EdgeColor','none')
        end
        hold off

        %labels
        title('Sorted by AvoidedEms - highest to lowest')
        xlim([0 CDFx(end,1)]);
        ylim([0 CDFy(end,1)]);
        xlabel('Production - kcal');
        ylabel('Annual emissions avoided - GtCO2');
    

    %Panel B - spared land

        % make CDF sorted by highest to lowest SparedLand
        Sorted = sortrows(Results,3,'descend');
        xvar = 2; %Production
        yvar = 3; %Area
        for i = 1:size(Results,1)
            if i==1
                CDFx(i,1) = Sorted(i,xvar);
                CDFy(i,1) = Sorted(i,yvar);
            else
                CDFx(i,1) = Sorted(i,xvar) + CDFx(i-1,1);
                CDFy(i,1) = Sorted(i,yvar) + CDFy(i-1,1);
            end
        end

        subplot(2,1,2)
        edges = zeros(length(CDFx)+1,1);
        edges(2:end,1) = CDFx;
        vals = CDFy;
        center = (edges(1:end-1) + edges(2:end))/2;
        width = diff(edges);
        hold on
        for i=1:length(center)
            bar(center(i),vals(i),width(i),'FaceColor',Sorted(i,4:6),'EdgeColor','none')
        end
        hold off

        %labels
        title('Sorted by SparedLand - highest to lowest')
        xlim([0 CDFx(end,1)]);
        ylim([0 CDFy(end,1)]);
        xlabel('Production - kcal');
        ylabel('Land area - ha');

% for text
PalmEms = sum(sum(AvoidedEms(:,:,1)));
SoyEms = sum(sum(AvoidedEms(:,:,2)));
disp('share of palm oil emissions:')
PalmEms/(PalmEms+SoyEms)

PalmArea = sum(sum(Area(:,:,1)));
SoyArea = sum(sum(Area(:,:,2)));
disp('share of palm oil area:')
PalmArea/(PalmArea+SoyArea)


    if PrintFigs == 1
       print(gcf,'-depsc','-painters', '../Plots/Figure3[CDFs]');
    end

end %if Fig3



