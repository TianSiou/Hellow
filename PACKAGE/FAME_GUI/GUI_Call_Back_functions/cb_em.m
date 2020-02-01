function Hline = cb_em(h,action)
global ith_wave_vector
global Pgraph
global hh
global hline 
Pgraph.wave_vec_array
Pgraph
% ith_wave_vector

if nargin < 2
    action='start';
end
    switch action
        case 'start'
            if isempty(hline) || ~ishandle(hline) 
                axes(h.ax_band);
                xL = xlim; yL = ylim;
                hold on
                hline = plot(mean(xL)*ones(2,1),yL,'r');
                xL
                ith_wave_vector = round(mean(xL))
                hold off
            end
            set(hline,'visible','on')
            set(gcf, 'WindowButtonDownFcn', 'cb_em([],''down'')');
            Hline = hline;
            hh    = h
        case 'down'
            currPt = get(gca, 'CurrentPoint');

            [~,ith_wave_vector] = min( abs(Pgraph.xdata - currPt(1,1)) )
hline
            temp = (Pgraph.wave_vec_array(:,ith_wave_vector))';
            set(hh.log_text(2),'string',['wave vector=[' num2str(temp) ']']);
            set(hline(1),'xdata',currPt(1,1)*ones(2,1));
            set(gcf, 'WindowButtonMotionFcn', 'cb_em([],''move'')');
            set(gcf, 'WindowButtonUpFcn', 'cb_em([],''up'')');
        case 'move'
            currPt = get(gca, 'CurrentPoint');
           % ith_wave_vector=round(currPt(1,1))+1;
            [~,ith_wave_vector]=min(abs(Pgraph.xdata-currPt(1,1)))
           % set(hh.log_text(2),'string',['wave vector=[' num2str(Pgraph.wave_vec_array(ith_wave_vector,:)) ']']);
            temp = (Pgraph.wave_vec_array(:,ith_wave_vector))';
            set(hh.log_text(2),'string',['wave vector=[' num2str(temp) ']']);
          
            set(hline,'xdata',currPt(1,1)*ones(2,1));
            
        case 'up'
            currPt = get(gca, 'CurrentPoint');
           % ith_wave_vector=round(currPt(1,1))+1;
            [~,ith_wave_vector]=min(abs(Pgraph.xdata-currPt(1,1)))
            temp = (Pgraph.wave_vec_array(:,ith_wave_vector))';
            set(hh.log_text(2),'string',['wave vector=[' num2str(temp) ']']);
          
           % set(hline,'xdata',(ith_wave_vector-1)*ones(2,1));
            set(hline,'xdata',Pgraph.xdata(ith_wave_vector)*ones(2,1));
           % xx
            set(gcf, 'WindowButtonMotionFcn', '');
            set(gcf, 'WindowButtonUpFcn', '');        
    end
%     ith_wave_vector
end