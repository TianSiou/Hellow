function Call_Back_save(hobj,event,h)
    [fn,path] = uiputfile('*.mat;*.bmp');
    if fn
        if strcmp('bmp',fn(end-2:end))
            set(h.fig,'PaperOrientation','portrait','PaperPosition',[0.5 0.5 13 9])
            f  = getframe(h.fig);
            im = frame2im(f);
            imwrite(im,[path fn(1:end-4) '.bmp'],'bmp');
        elseif strcmp('mat',fn(end-2:end))
%             try
                load savedata
                save([path fn],'Popt','Pgraph','Pmesh','Pgrid','Pnum','Pvec','Pmaterial','Plattice','par');
%             end
        end
    end
end