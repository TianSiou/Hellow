function FAME_Plot_Band_Structure( path_string, part_num, Freq_array )
figure(1)
if iscell(path_string) == 0
    flag = 0;
    n_new = 1;
    for i = 1:length(path_string)-1
        if strcmp(path_string(i),'|') == 1
            path_string_new{n_new} = path_string(i-1:i+1);
            n_new = n_new + 1;
            flag  = 1;
        elseif strcmp(path_string(i+1),'|') == 1
        elseif flag == 1
            flag = 0;
        else
            path_string_new{n_new} = path_string(i);
            n_new = n_new + 1;
        end
    end
    path_string_new{n_new} = path_string(i+1);
    path_string = path_string_new;
end

% Freq_array = sort(Freq_array,'ascend');
%% Plotting the band structure
% Starting to plot the data 
plot_num     = size(Freq_array,1);
node_num     = length(path_string);
segment_pt   = 1;
segment_len  = 0;
for i = 1:length(path_string)
    if length(path_string{i}) == 3 ||  i == length(path_string)
        segment_pt  = [segment_pt,i];
        segment_len = [ segment_len, (i - segment_pt(end-1))*(part_num-1) + 1 ];
    end
end
% xdata = [];
% for j = 1:length(segment_len)-1
for j = 1:length(segment_len)-1
%     temp = [];
    for plot_idx = 1 : plot_num   
       x_temp = sum(segment_len(1:j))-j+1:sum(segment_len(1:j))+segment_len(j+1)-j;
%        temp = [temp, x_temp];
       y_temp = Freq_array( plot_idx , sum(segment_len(1:j))+1:sum(segment_len(1:j+1)) );
       plot( x_temp.' , y_temp,'color', [ 0.1 0.05 0.7 ],'linewidth', 1.5 );
%         plot( x_temp.' , y_temp,'color', [ 1 0 0 ],'linewidth', 1.5 );
       hold on
    end
%     xdata(j,:)  = temp
end
% Setting the labels, axis and grid line on the band structure graph 
title ('Band Structure','fontsize',14 );
xlabel('Wave Vector'   ,'fontsize',11 ); 
ylabel('Frequency'     ,'fontsize',11 ); 
grid on;
axis tight;
label_symbol_array = cell ( 1 , node_num ); 
label_number_array = zeros( 1 , node_num );
label_symbol_array(1) = {path_string{1}};
if node_num > 1
    for path_idx = 2 : node_num 
        label_number_array(path_idx) = label_number_array(path_idx-1) + part_num-1;
        label_symbol_array(path_idx) = {path_string{path_idx}} ;
        if length(path_string{path_idx}) == 3
            plot([label_number_array(path_idx),label_number_array(path_idx)],[0,max(max(Freq_array))/(2*pi)],'k--')
        end
    end
end
set(gca,'FontSize'  ,                         12);
set(gca,'XTick'     ,       label_number_array  );
set(gca,'XTickLabel',       label_symbol_array );    