function Call_Back_clear(hobj,event,h)
% 清除所有設定, 回到初始設定
    button = questdlg('Are you sure to reset all the parameters?','','Yes','Save as...','Cancel','Yes');
    
    switch button
        case 'Yes'
            FAME_GUI;
        case 'Cancel'
            
        case 'Save as...'
            Call_Back_save([],[],h);
    end
end