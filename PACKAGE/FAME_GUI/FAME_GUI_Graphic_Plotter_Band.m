function FAME_GUI_Graphic_Plotter_Band(  eigvalue_array, path_string, part_num, xdata,xtick,xticklabal)
%% Plotting the band structure
figure(1)
% Starting to plot the data 
wv  =  xdata;

for plot_idx = 1 : size(eigvalue_array,1)
    plot( wv , eigvalue_array( plot_idx , : ),'bo-','linewidth', 1.5 )   
    hold on
end  
% Setting the labels, axis and grid line on the band structure graph 
title ('Band Structure','fontsize',18 );
xlabel('Wave Vector'   ,'fontsize',18 ); 
ylabel('Frequency'     ,'fontsize',18 ); 
grid off;
axis tight;

% label_number_array = zeros( 1 , node_num );

set(gca,'FontSize'  ,                  16);
set(gca,'XTick'     ,       xtick );
set(gca,'XTickLabel', upper(xticklabal) );

end


