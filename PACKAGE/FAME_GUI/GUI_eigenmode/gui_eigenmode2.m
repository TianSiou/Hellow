function gui_eigenmode2(action,eigvector,eigvalue,SN)
    persistent haxes hgca hobj htext hslider hpush idx_gca htimer hpopup hfig
    persistent xyz xyz_old b1 b2 flag x y z himg idx_xyz cx cy cz flag_check hcheck hslice hcones tt flag_doubleclick minuvw maxuvw range range_time rrange num_eig neig
    persistent hquiver scalefactor eigenvector eigenvalue
    persistent u v w u0 v0 w0 Time_harmonic_component sn ss
    global ln Pgraph
%%%%%%%%
%ss=2;
    warning off
    if nargin==4
        sn = SN;%[16,16,16,3];
    end
    range=[0 1
           0 1
           0 1];
    if nargin > 2
        eigenvector = eigvector;
        eigenvalue  = eigvalue;
    end
    if isempty(ss)
        ss=2;
    end
    if isempty(flag_doubleclick)
        flag_doubleclick=1;
    end
    if nargin<1
        action='start';
    end
    time_interval=0.05;
    idouble_pi=2*pi*1i;

    rrange=range';
    aaxis=[ range(2,:)  range(1,:) 
            range(3,:)  range(1,:) 
            range(2,:)  range(3,:) ];
    idx=[1 2; 1 3; 3 2];
    cc=['brb';'kkr'];
    % range_time=[0 20];   
%%%%%%%%    
    switch action
        case 'start'
            num_eig = length(eigenvalue);
            neig    = 1;
            scalefactor = 1;         

            Time_harmonic_component = exp( -idouble_pi * eigenvalue(neig) * ( 0 : time_interval : (1 / eigenvalue(neig) ) ) );
            ln=length(Time_harmonic_component); 
            range_time = [0 ln-1]; 
            [x,y,z] = meshgrid(range(1,1):(range(1,2)-range(1,1))/(sn(1)-1):range(1,2),...
                               range(2,1):(range(2,2)-range(2,1))/(sn(2)-1):range(2,2),...
                               range(3,1):(range(3,2)-range(3,1))/(sn(3)-1):range(3,2));
            xmin = min(x(:));
            xmax = max(x(:));
            ymin = min(y(:));
            ymax = max(y(:));
            zmin = min(z(:));
            zmax = max(z(:));

            uvw = reshape(eigenvector(:,neig),sn);%*Time_harmonic_component(1));
% maxuvw=max(uvw(:));
% minuvw=min(uvw(:));

            xrange = linspace(range(1,1),range(1,2),8);
            yrange = linspace(range(2,1),range(2,2),8);
            zrange = linspace(range(3,1),range(3,2),8);
            [cx cy cz] = meshgrid(xrange,yrange,zrange);
%             v = ux.*exp(-x.^2-y.^2-z.^2);
%             u = v;...x.*exp(+x.^2+y.^2+z.^2);
%             w=v;
            u0 = uvw(:,:,:,1);
            v0 = uvw(:,:,:,2);
            w0 = uvw(:,:,:,3);
            u  = 5*real(u0*Time_harmonic_component(1));
            v  = 5*real(v0*Time_harmonic_component(1));
            w  = 5*real(w0*Time_harmonic_component(1));
            maxuvw = max(abs([u(:);v(:);w(:)]));          
            minuvw = -maxuvw;             
            load bb
            flag = 1;
            idx = [1 2; 1 3; 3 2];
            xyz = [0.1 .5 0.5];xyz_old=xyz;
                    %idx_xyz=ceil((xyz-[0,0,0])*15);
                    %idx_xyz=floor((xyz-range(:,1)')./(range(:,2)'./(sn(1:3)-1)))+1;
            idx_xyz = ceil((xyz-range(:,1)')./(range(:,2)'./(sn(1:3)-1)));
            
            % 開啟視窗
            bcolor = [1 1 1 ]/1.06;
            if isempty(hfig) || ~ishandle(hfig)
                hfig = figure( 'Name',        ['Eigenmode: wave vector=[' num2str(Pgraph.wave_vec_array(:,Pgraph.ith_wave_vector)') ']'] ,...
                               'NumberTitle', 'off',...
                               'DeleteFcn',   sprintf('%s %s', mfilename, 'close'));
                clf(hfig);
            else
                set(hfig,'Name','Eigenmode','NumberTitle','off','DeleteFcn',sprintf('%s %s', mfilename, 'close'));
                clf(hfig);
            end
            set(hfig, 'position', [0 0 1200 800],...
                      'color',    bcolor,...
                      'menubar', 'none',...
                      'toolbar', 'none');
            % 給定各切面 axis      
%             pos = [.025 .6 .45 .38
%                     .51 .6 .45 .38                
%                    .025 .2 .45 .38
%                     .51 .2 .45 .38];
            pos = [.025 .6 .3 .3
                    .31 .6 .3 .3                
                   .025 .2 .3 .3
                    .31 .2 .3 .3];
            pos(:,3:4) = pos(:,3:4)-0.04;
            pos(:,1:2) = pos(:,1:2)+0.04;
            for i = 1:3
                haxes(i) = axes('parent',hfig,'position',pos(i,:));
                axes(haxes(i));
                       
                if i == 1 %xy
                    himg(i) = imagesc(y(:,1,1),squeeze(x(1,:,1)),squeeze(w(:,:,idx_xyz(3))),[minuvw maxuvw]);
                    hold on
                    hquiver(i) = quiver(haxes(i),y(1:ss:end,1,1),x(1,1:ss:end,1),scalefactor*squeeze(v(1:ss:end,1:ss:end,idx_xyz(3))),scalefactor*squeeze(u(1:ss:end,1:ss:end,idx_xyz(3))),0,'color',[1 1 1],'LineWidth',2);
                    hold off
                    xlabel('y','color','r','fontsize',12);ylabel('x','color','r','fontsize',12);

                elseif i == 3 %zy
                    himg(i) = imagesc(squeeze(y(:,1,1)),squeeze(z(1,1,:)),squeeze(u(idx_xyz(1),:,:))',[minuvw maxuvw]);  
                    hold on
                    hquiver(i)=quiver(squeeze(y(1:ss:end,1,1)),squeeze(z(1,1,1:ss:end)),scalefactor*squeeze(v(idx_xyz(1),1:ss:end,1:ss:end)),scalefactor*squeeze(w(idx_xyz(1),1:ss:end,1:ss:end)),0,'color',[1 1 1],'LineWidth',2);
                    hold off
                    xlabel('y','color','r','fontsize',12);ylabel('z','color','r','fontsize',12);
                else %xz
                    himg(i) = imagesc(squeeze(z(1,1,:)),squeeze(x(1,:,1)),squeeze(v(:,idx_xyz(2),:)),[minuvw maxuvw]);
                    hold on
                    hquiver(i)=quiver(squeeze(z(1,1,1:ss:end)),x(1,1:ss:end,1),scalefactor*squeeze(w(1:ss:end,idx_xyz(2),1:ss:end)),scalefactor*squeeze(u(1:ss:end,idx_xyz(2),1:ss:end)),0,'color',[1 1 1],'LineWidth',2);
                    hold off
                    xlabel('z','color','r','fontsize',12);ylabel('x','color','r','fontsize',12);
                end
                hold on 
                hobj(i,1)=plot(xyz(idx(i,2))*ones(1,2),aaxis(i,3:4),cc(1,i)); 
                hobj(i,2)=plot(aaxis(i,1:2),xyz(idx(i,1))*ones(1,2),cc(2,i)); 
 % hgtext(i)=text(xyz(idx(i,2)),xyz(idx(i,1)),num2str(xyz));
                hold off
%axis(aaxis(i,:));  %axis ij
                set( haxes(i),'XColor',[1 1 1],'YColor',[1 1 1]); 
                axis equal;axis(aaxis(i,:));
                drawnow                                     
            end
            hgca    = haxes(1);
            idx_gca = 1;
            set(haxes(1),'YDir','reverse')
            set(haxes(2),'YDir','reverse')
            set(haxes(3),'YDir','normal')
            i = 4;
            haxes(i) = axes('parent',hfig,'position',pos(i,:));

            % 顯示時間用 slider
            hslider = uicontrol(hfig, 'style',      'slider',...
                                      'units',      'normalized',...
                                      'position',   [.12 .1 .44 .05],...
                                      'max',        range_time(2),...
                                      'min',        range_time(1),...
                                      'SliderStep', [0.01 0.01],...
                                      'value',      0, ...
                                      'callback',   sprintf('%s %s', mfilename, 'cb_slider'));              
            hslider(2) = uicontrol(hfig, 'style',      'slider',...
                                         'units',      'normalized',...
                                         'position',   [.60 .83 .38 .05],...
                                         'max',        1,...
                                         'min',        0,...
                                         'SliderStep', [0.05 0.05],...
                                         'value',      1, ...
                                         'callback',   sprintf('%s %s', mfilename, 'cb_slider2'));
            hslider_text(2) = uicontrol(hfig, 'style',           'text',...
                                              'units',           'normalized',...
                                              'position',        get(hslider(2),'position') + [0 .05 -.2 -.01],...
                                              'string',          'transparency',...
                                              'fontsize',        16,...
                                              'HorizontalAlignment', 'left',...
                                              'ForegroundColor', [0 0 1],...
                                              'backgroundcolor', bcolor);
            hslider(3) = uicontrol(hfig, 'style',      'slider',...
                                         'units',      'normalized',...
                                         'position',   [.6 .73 .38 .05],...
                                         'max',        40,...
                                         'min',        0,...
                                         'SliderStep', [0.05 0.05],...
                                         'value',      1, ...
                                         'callback',   sprintf('%s %s', mfilename, 'cb_slider3'));
            hslider_text(3) = uicontrol(hfig, 'style',           'text',...
                                              'units',           'normalized',...
                                              'position',        get(hslider(3),'position') + [0 .05 -.2 -.01],...
                                              'string',          'arrow size',...
                                              'fontsize',        16,...
                                              'HorizontalAlignment', 'left',...
                                              'ForegroundColor', [0 0 1],...
                                              'backgroundcolor', bcolor);
            % 播放用的 push
            hpush = uicontrol(hfig, 'style',    'push',...
                                    'units',    'normalized',...
                                    'position', [.05 .1 .06 .05],...
                                    'cdata',    b2,...
                                    'callback', sprintf('%s %s', mfilename, 'cb_push'));
            htext = uicontrol(hfig, 'style',               'text',...
                                    'string',              [],...
                                    'units',               'normalized',...
                                    'fontsize',            16,...
                                    'HorizontalAlignment', 'left',...
                                    'foregroundcolor',     [0 0 0],...
                                    'backgroundcolor',     bcolor,...
                                    'position',            [.6 .1 .4 .04]);
            axes(haxes(4));hold on;
            
            set( haxes(4),'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]); 
            hslice(3) = slice(x,y,z,w,[],[],xyz(3));hold on
            hslice(1) = slice(x,y,z,u,xyz(1),[],[]);
            hslice(2) = slice(x,y,z,v,[],xyz(2),[]);
            %set(hslice,'EdgeAlpha',0,'FaceAlpha',1);
            set(hslice, 'EdgeAlpha',0,...get(hslider(2),'value'),
                        'FaceAlpha',get(hslider(2),'value'));
            set(haxes(4),'CLim',[minuvw maxuvw]);
            hcones = coneplot(x,y,z,scalefactor*u,scalefactor*v,scalefactor*w,cx,cy,cz,0);
            set(hcones,'FaceColor',[238 201 0]/255,'EdgeColor','none')
            hold off
            %xlabel('x');ylabel('y');zlabel('z');
            view(52,28)    
            axis equal;axis(rrange(:));
            light;light;lighting phong
            drawnow
            %  xlabel('x','color','r','fontsize',12);ylabel('y','color','r','fontsize',12);zlabel('z','color','r','fontsize',12);

            set(gcf, 'WindowButtonDownFcn', sprintf('%s %s', mfilename, 'down'));
            htimer = timer('timerfcn',{@slidershow,hslider(1)},'ExecutionMode','fixedSpacing','period',.001);
            
            hpopup_text(1) = uicontrol(hfig,  'style',           'text',...
                                              'units',           'normalized',...
                                              'position',        [.6 .63 .14 .04],...
                                              'string',          'color map: ',...
                                              'fontsize',        16,...
                                              'HorizontalAlignment', 'left',...
                                              'ForegroundColor', [0 0 1],...
                                              'backgroundcolor', bcolor);
            hpopup(1) = uicontrol(hfig, 'style',           'popup',...
                                        'string',          'Jet|HSV|Hot|Cool|Gray|Copper',...
                                        'units',           'normalized',...
                                        'fontsize',        16,...
                                        'foregroundcolor', [0 0 1],...
                                        'backgroundcolor', [1 1 1],...
                                        'position',        get(hpopup_text(1), 'position') + [.15 0 0 0],...
                                        'callback',        sprintf('%s %s', mfilename, 'cb_hpopup'));
              
            cc=[];for i=1:num_eig; cc=[cc int2str(i) '|'];end;cc(end)=[];
            hpopup_text(2) = uicontrol(hfig,  'style',           'text',...
                                              'units',           'normalized',...
                                              'position',        [.6 .53 .14 .04],...
                                              'string',          '# of eigenvector: ',...
                                              'fontsize',        16,...
                                              'HorizontalAlignment', 'left',...
                                              'ForegroundColor', [0 0 1],...
                                              'backgroundcolor', bcolor);
            hpopup(2) = uicontrol(hfig, 'style',           'popup',...
                                        'string',          cc,...
                                        'units',           'normalized',...
                                        'fontsize',        16,...
                                        'foregroundcolor', [0 0 1],...
                                        'backgroundcolor', [1 1 1],...
                                        'position',        get(hpopup_text(2), 'position') + [.15 0 0 0],...
                                        'callback',        sprintf('%s %s', mfilename, 'cb_hpopup2')); 
                                    
            cc=[];for i=1:4; cc=[cc int2str(2^(i-1)) '|'];end;cc(end)=[];
            hpopup_text(3) = uicontrol(hfig,  'style',           'text',...
                                              'units',           'normalized',...
                                              'position',        [.6 .43 .14 .04],...
                                              'string',          'sparsity of arrow: ',...
                                              'fontsize',        16,...
                                              'HorizontalAlignment', 'left',...
                                              'ForegroundColor', [0 0 1],...
                                              'backgroundcolor', bcolor);
            hpopup(3) = uicontrol(hfig, 'style',           'popup',...
                                        'string',          cc,...
                                        'units',           'normalized',...
                                        'fontsize',        16,...
                                        'foregroundcolor', [0 0 1],...
                                        'backgroundcolor', [1 1 1],...
                                        'position',        get(hpopup_text(3), 'position') + [.15 0 0 0],...
                                        'value',           log2(ss)+1,...
                                        'callback',        sprintf('%s %s', mfilename, 'cb_hpopup3'));                                     
            cc={'xy','yz','xz','cone','cross'};                         
            flag_check=[1 1 1 1 1];
            for i = 1:5
                hcheck(i) = uicontrol(hfig, 'style',           'check',...
                                            'string',          cc{i},...
                                            'units',           'normalized',...
                                            'fontsize',        16,...
                                            'foregroundcolor', [0 0 1],...
                                            'backgroundcolor', bcolor,...
                                            'position',        [.025+i*.09 .95 .1 .04],...
                                            'value',           1,...
                                            'callback',        sprintf('%s %s', mfilename, 'cb_check'));
                                        
            end
%        set(Pgraph.par.hsc,'parent',haxes(4))
        case 'down'
           % currPt1 = get(gca, 'CurrentPoint');
            hgca=gca;
            tt=tic;
            for i=1:4
                if hgca==haxes(i)
                    idx_gca=i;
                    break;
                end                  
            end
            set(gcf, 'WindowButtonMotionFcn', sprintf('%s %s', mfilename, 'move'));
            set(gcf, 'WindowButtonUpFcn', sprintf('%s %s', mfilename, 'up'));
        case 'move'  
            eval(sprintf('%s %s', mfilename, 'showimage'));
        case 'up'
            eval(sprintf('%s %s', mfilename, 'oneclick'));
            eval(sprintf('%s %s', mfilename, 'showimage'));
            set(gcf, 'WindowButtonMotionFcn', '');
            set(gcf, 'WindowButtonUpFcn', '');
        case 'cb_slider'  
            set(htext,'string',['[x y z]=[' num2str(xyz) '], time=' num2str(get( hslider(1),'value'))]);
            ln = length(Time_harmonic_component);
            t  = fix(get(hslider(1),'value'));
             %set(h,'value',mod(t+1,ln));
             %gui_eigenmode2('cb_slider');
             %     u=u0*real(exp(j*t*pi));
             %    v=v0*real(exp(j*t*pi));
             %    w=w0*real(exp(j*t*pi));
            u = 5*real(u0*Time_harmonic_component(mod(t,ln)+1));
            v = 5*real(v0*Time_harmonic_component(mod(t,ln)+1));
            w = 5*real(w0*Time_harmonic_component(mod(t,ln)+1));
            gui_eigenmode2('showimage');drawnow
        case 'cb_slider2'  
            set(hslice, 'EdgeAlpha', 0,...get(hslider(2),'value'),
                        'FaceAlpha', get(hslider(2),'value'));
        case 'cb_slider3'
            scalefactor = get(hslider(3),'value');
            eval(sprintf('%s %s', mfilename, 'showimage'));
        case 'cb_push'  
            if flag
                flag=0;
                set(hpush,'cdata',b1);
                start(htimer);
            else
                flag=1;
                set(hpush,'cdata',b2);
                stop(htimer);
            end
        case 'cb_check'
            for i=1:5
                flag_check(i)=get(hcheck(i),'value');
            end
            cc={'off','on'};oidx=[2 3 1];
            for i=1:3
                if ishandle(hslice(i))
                    set(hslice(i),'visible',cc{flag_check(oidx(i))+1});
                end
            end
            if ishandle(hcones)
                set(hcones,'FaceColor',[238 201 0]/255,'EdgeColor','none','visible',cc{flag_check(4)+1});
            end    
            set(hobj,'visible',cc{flag_check(5)+1});
        case 'cb_hpopup'
            cmap = get(hpopup(1),'string');
            v    = get(hpopup(1),'value');
            eval(['colormap ' cmap(v,:)]);
        case 'cb_hpopup2'
            neig = get(hpopup(2),'value');
            %%
            Time_harmonic_component = exp( -idouble_pi * eigenvalue(neig) * ( 0 : time_interval : (1 / eigenvalue(neig) ) ) );
            ln=length(Time_harmonic_component); 
            range_time=[0 ln-1]; 
            
            uvw = reshape(eigenvector(:,neig),sn);%*Time_harmonic_component(1));
            u0  = uvw(:,:,:,1);
            v0  = uvw(:,:,:,2);
            w0  = uvw(:,:,:,3);
            maxuvw = max([u(:);v(:);w(:)]);          
            minuvw = min([u(:);v(:);w(:)]); 
            gui_eigenmode2('cb_slider');
        case 'cb_hpopup3'
            ss = 2^(get(hpopup(3),'value')-1);
            gui_eigenmode2('start',eigenvector,eigenvalue,sn);
        case 'showimage'
            currPt = get(hgca, 'CurrentPoint');
        %  if idx_gca <4
            if idx_gca <4 && currPt(1,1)>=aaxis(idx_gca,1) && currPt(1,1) <=aaxis(idx_gca,2) &&  currPt(1,2)>=aaxis(idx_gca,3) && currPt(1,2) <=aaxis(idx_gca,4)
                xyz(  idx(idx_gca,:))=  currPt(1,[2 1]);
                xyz_old=xyz;
            else
                xyz=xyz_old;
            end
           %idx_xyz=floor((xyz-range(:,1)')./(range(:,2)'./(sn(1:3)-1)))+1;
            idx_xyz=ceil((xyz-range(:,1)')./(range(:,2)'./(sn(1:3)-1)));
            for i=1:3
                if i==1
                    set(himg(i),'cdata',squeeze(w(:,:,idx_xyz(3)))');
                    set(hquiver(i),'udata',scalefactor*squeeze(v(1:ss:end,1:ss:end,idx_xyz(3)))','vdata',scalefactor*squeeze(u(1:ss:end,1:ss:end,idx_xyz(3))'),'LineWidth',2);   
                elseif i==3
                    set(himg(i),'cdata',squeeze(u(idx_xyz(1),:,:))');
                    set(hquiver(i),'udata',scalefactor*squeeze(v(idx_xyz(1),1:ss:end,1:ss:end))','vdata',scalefactor*squeeze(w(idx_xyz(1),1:ss:end,1:ss:end))','LineWidth',2);
                else
                    set(himg(i),'cdata',squeeze(v(:,idx_xyz(2),:))); 
                    set(hquiver(i),'udata',scalefactor*squeeze(w(1:ss:end,idx_xyz(2),1:ss:end)),'vdata',scalefactor*squeeze(u(1:ss:end,idx_xyz(2),1:ss:end)),'LineWidth',2);
                end
                set(  hobj(i,1),'xdata',xyz(idx(i,2))*ones(1,2),'ydata', aaxis(i,3:4));          
                set(  hobj(i,2),'xdata',aaxis(i,1:2),'ydata', xyz(idx(i,1))*ones(1,2));         
                set(haxes(i),'CLim',[minuvw maxuvw]);
            end
            set(htext,'string',['[x y z]=[' num2str(xyz) '], time=' num2str(get( hslider(1),'value'))]);drawnow
            axes(haxes(4));cla; hold on; 
            cc = {'off','on'};oidx=[2 3 1];
            if flag_check(oidx(3))
                hslice(3)=slice(x,y,z,w,[],[],xyz(3),'nearest');
            end
            if flag_check(oidx(1))
                hslice(1)=slice(x,y,z,u,xyz(1),[],[],'nearest');
            end
            if flag_check(oidx(2))
                hslice(2)=slice(x,y,z,v,[],xyz(2),[],'nearest');
            end
            for i = 1:3
                if ishandle(hslice(i))
                    set(hslice(i),'EdgeAlpha', 0,...
                                  'FaceAlpha', get(hslider(2),'value'));
                end
            end
            cc = {'off','on'};oidx=[2 3 1];
            if flag_check(4)
                hcones = coneplot(x,y,z,scalefactor*u,scalefactor*v,scalefactor*w,cx,cy,cz,0);
                set(hcones,'FaceColor',[238 201 0]/255,'EdgeColor','none','visible',cc{flag_check(4)+1});
            end
            set( haxes(4),'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
            hold off
            view(52,28);light;light;lighting phong
            axis equal;axis(rrange(:));set(haxes(4),'CLim',[minuvw maxuvw]);
%            xlabel('x','color','r','fontsize',12);ylabel('y','color','r','fontsize',12);zlabel('z','color','r','fontsize',12);
            drawnow   
        case 'oneclick'
%             pos = [.025 .6 .45 .38
%                     .51 .6 .45 .38                
%                    .025 .2 .45 .38
%                     .51 .2 .45 .38];
            pos = [.025 .6 .3 .3
                    .31 .6 .3 .3                
                   .025 .2 .3 .3
                    .31 .2 .3 .3];
            if idx_gca == 4
                if toc(tt) < .5 
                    if flag_doubleclick 
                        for i = 1:3
                            set(haxes(i),'position',[0 0 0.01 0.01],'visible','off');
                        end
                        for i = 1:5
                            t = i/5;
                            set(haxes(4),'position',(1-t)*pos(4,:)+t*[.001 .2 .7 .7]);
                            drawnow;% pause(.001)
                        end
                        %rotate3d on
                        flag_doubleclick=0;
                    else
                        for i = 1:3
                            set(haxes(i),'position',pos(i,:));
                        end
                        for i = 1:5
                            t = i/5;
                            set(haxes(4),'position',(t)*pos(4,:)+(1-t)*[.001 .2 .7 .7]);
                            drawnow
                        end
                        %rotate3d off
                        flag_doubleclick=1;
                    end
                end
            end
        case 'close'
            if strcmp(get(htimer,'running'),'on')
                stop(htimer);
                delete(htimer);
            end
            close(hfig)
    end
end

function slidershow(obj,event,h)
    global ln %Time_harmonic_component

    t=fix(get(h,'value'));
    set(h,'value',mod(t+1,ln));
    gui_eigenmode2('cb_slider');
end






