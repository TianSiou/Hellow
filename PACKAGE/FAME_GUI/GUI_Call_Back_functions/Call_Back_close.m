function Call_Back_close(hobj,event,h)
% �M���Ҧ��]�w, �^���l�]�w
    button = questdlg('Are you sure to exit this program?','','Yes','Save as...','Cancel','Yes');
    
    switch button
        case 'Yes'
            close all
        case 'Cancel'
            
        case 'Save as...'
            Call_Back_save([],[],h);
    end
      
end