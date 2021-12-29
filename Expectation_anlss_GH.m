% ========================== Opening / closing functions =================

function varargout = Expectation_anlss_GH(varargin)

    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Expectation_anlss_GH_OpeningFcn, ...
                       'gui_OutputFcn',  @Expectation_anlss_GH_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    
end

function Expectation_anlss_GH_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    
    handles.main_folder = 'G:\Expectation' ;
    
    set(handles.figure1,'windowscrollWheelFcn',@scroll_func)
        
    handles.phys_mice = [1:3 5] ;
    
    handles.prob_colors = {[0.474 0.623 0.796],[0.5 0.5 0.5],[0.976 0.4 0.369]} ;
    handles.exp_colors = {[0 0.659 0.643],[0.5 0.5 0.5],[0.651 0.549 0.89]} ;
                
    handles.data_holder = {} ;

    handles.load_data_button.Enable = 'on' ;
    
    guidata(hObject,handles) ;
    
end

function varargout = Expectation_anlss_GH_OutputFcn(hObject, eventdata, handles) 

    varargout{1} = handles.output;
    
end

% ========================== Main callbacks & functions ==================

function all_mice_anlss_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    if handles.save_anlss_ax.UserData
        figure() ;
        Ax = gca ;
    else
        Ax = handles.anlss_ax ;
    end

    handles.all_mice_anlss_menu.UserData = 1 ;
    handles.one_mouse_anlss_menu.UserData = 0 ;
    menu = handles.all_mice_anlss_menu ;
    cla(Ax , 'reset')
    axis(Ax , 'xy')
    axis(Ax , 'auto')
    hold(Ax , 'on')
    set(Ax,'Visible','on','Box','off')
    set(Ax, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')

    switch menu.String{menu.Value}
        case 'Bias progress'
            plot_odor_behav_prog(Ax,hObject,'bias')
        case 'FA rate progress'
            plot_odor_behav_prog(Ax,hObject,'FA rate')
        case 'Hit rate progress'
            plot_odor_behav_prog(Ax,hObject,'hit rate')
        case "d' expectation progress"
            plot_odor_behav_prog(Ax,hObject,"d' expectation")
        case "d' progress"
            plot_gen_behav_prog(Ax,hObject,"d'")
        case 'Selected days bias'
            plot_select_days_behav(Ax,hObject,'bias')
        case "Selected days d' expectation"
            plot_select_days_behav(Ax,hObject,"d' expectation")
        case 'Selected days HR'
            plot_select_days_behav(Ax,hObject,'HR')
        case 'Selected days FAR'
            plot_select_days_behav(Ax,hObject,'FAR')
        case 'sound vs. odor'
            plot_odor_sound_behav_prog(Ax,hObject)
        case 'Expectation DI (bias)'
            all_auditory_neurons(Ax,hObject,'expectation DI','bias','violin')
        case 'DI progress'
            plot_DI_progress(Ax,hObject)
        case 'DI diff. progress'
            plot_DI_diff_progress(Ax,hObject)
        case "dDI to d'"
            plot_dDI_to_dprime(Ax,hObject)
        case 'Population mean (bias)'
            plot_all_population_means(hObject,'bias')

    end
        

end

function one_mouse_anlss_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    if handles.save_anlss_ax.UserData
        figure() ;
        Ax = gca ;
    else
        Ax = handles.anlss_ax ;
    end
    
    handles.all_mice_anlss_menu.UserData = 0 ;
    handles.one_mouse_anlss_menu.UserData = 1 ;

    menu = handles.one_mouse_anlss_menu ;
    cla(Ax,'reset')
    axis(Ax , 'xy')
    axis(Ax , 'auto')
    hold(Ax , 'on')
    set(Ax,'Visible','on','Box','off')
    set(Ax, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')

    switch menu.String{menu.Value}
        case 'Neuronal field'
            display_field(Ax,hObject)
        case 'Trial counts'
            plot_trial_counts(Ax,hObject)
        case 'Lick rate'
            plot_lick_rate(Ax,hObject)
        case "d'"
            plot_dprime(Ax,hObject) 
        case 'Raw behavior'
            plot_raw_behav(hObject) 
        case 'Single session behavior'
            plot_single_sess_behav_all(Ax,hObject)
        case 'Single session HR'
            plot_single_sess_behav_parts(Ax,hObject,'HR')
        case 'Single session FAR'
            plot_single_sess_behav_parts(Ax,hObject,'FAR')
        case 'Single session bias'
            plot_single_sess_behav_parts(Ax,hObject,'bias')
        case "Single session d'"
            plot_single_sess_behav_parts(Ax,hObject,"d' expectation")
        case 'All sessions bias'
            plot_all_sess_behav(Ax,hObject,'bias')
        case 'All sessions FA rate'
            plot_all_sess_behav(Ax,hObject,'FA rate')                    
        case 'All sessions hit rate'
            plot_all_sess_behav(Ax,hObject,'hit rate')
        case 'Plot all neurons'
            plot_all_neurons(hObject)

    end
    
end

function update_anlss_menu(hObject)

    handles = guidata(hObject) ;

    handles.all_mice_anlss_menu.String = {'Bias progress',...
        'FA rate progress','Hit rate progress',...
        "d' expectation progress",'sound vs. odor'...
        'Selected days bias',"Selected days d' expectation",...
        'Selected days HR','Selected days FAR',...
        'Expectation DI (bias)','DI progress','DI diff. progress',...
        "dDI to d'",'Population mean (bias)'} ;

    handles.one_mouse_anlss_menu.String = {'Neuronal field' , ...
        'Single session behavior','Single session HR',...
        'Single session FAR','Single session bias',"Single session d'",'All sessions bias',...
        'All sessions FA rate','All sessions hit rate',...
        'Raw behavior','Trial sequence','Plot all neurons'} ;
    
    handles.all_mice_anlss_menu.Value = 1 ;
    handles.one_mouse_anlss_menu.Value = 1 ;
        
end

function neuron_filter_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;

    mouse = handles.mouse_menu.Value ;
    sess = handles.session_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;

    handles.neuron_menu.String = {} ;
    for neur = 1  : h.n_neurons
        neuron_crit = 0 ; % Filter criteria
        n = h.Neurons(neur) ;
        switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
            case 'All neurons'
                neuron_crit = 1 ; 
            case 'Responsive neurons'
                if ~n.invalid && n.responsive
                    neuron_crit = 1 ;
                end
            case 'Go neurons'
                if ~n.invalid && n.responsive && n.pref_sound==1
                    neuron_crit = 1 ;
                end
            case 'NoGo neurons'
                if ~n.invalid && n.responsive && n.pref_sound==0
                    neuron_crit = 1 ;
                end

        end
        if neuron_crit
            handles.neuron_menu.String{end+1} = ['Neuron ' num2str(neur)] ;
        end
    end
    if isempty(handles.neuron_menu.String)
        handles.neuron_menu.String = {'No neurons'} ;
    end
    
    if ~strcmp(handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value},'All neurons')
        handles.neuron_menu.Value = 1 ;
    end
    
    neuron_menu_Callback(hObject,0,0)

end

function update_neuron_filter_menu(hObject)

    handles = guidata(hObject) ;

    handles.neuron_filter_menu.String = {'All neurons','Responsive neurons',...
        'Go neurons','NoGo neurons'} ;
    handles.neuron_filter_menu.Value = 1 ;
    
end

% ========================== Single mouse analysis =======================

function display_field(Ax,hObject)

    handles = guidata(hObject) ;
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    legend(Ax,'hide')
    
    set(Ax,'YLimMode','auto','XLimMode','auto')
    xlabel(Ax,'')
    ylabel(Ax,'')

    col_hmask=sum(h.mask,3); % Collapsed hmask
    ilu_hmask=[zeros(size(col_hmask,1),1) diff(col_hmask,2,2) zeros(size(col_hmask,1),1)]+[zeros(1,size(col_hmask,2)) ; diff(col_hmask,2,1) ; zeros(1,size(col_hmask,2))]; % illustrated hmask

    field_with_mask = h.field;
    field_with_mask(:,:,3)=ilu_hmask;
    imshow(field_with_mask,'Parent',Ax)
    axis(Ax,'square')

    for i = 1:size(h.mask,3)
        temp = regionprops(h.mask(:,:,i),'centroid') ;
        text(Ax,temp.Centroid(1),temp.Centroid(2),num2str(i),'fontsize',14,'color','blue')
    end
    
end

function plot_accuracy(Ax,hObject)

    handles = guidata(hObject) ;
    legend(handles.anlss_ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    hit = [] ;
    miss = [] ;
    CR = [] ;
    FA = [] ;
    accuracy = [] ;
    
    for i = h.ft:10:h.lt
        hit(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==1) ;
        CR(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==2) ;
        miss(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==3) ;
        FA(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==4) ;
        n = length(hit) ;
        m = max([1 n-9]) ;
        accuracy(end+1) = 100*(sum(hit(m:n)) + sum(CR(m:n)))/...
                          (sum(hit(m:n)) + sum(CR(m:n)) + sum(miss(m:n)) + sum(FA(m:n))) ;
    end 
    
    hold(Ax , 'on')
    plot(Ax , 1:length(accuracy) , accuracy , 'k') ;
    plot(Ax , 1:length(accuracy) , accuracy , '*k')
    ylim(Ax , [50 100]) 
    
    xlabel(Ax,'Trials/10','fontsize',16)
    ylabel(Ax,'Accuracy (%)','fontsize',16)

end

function plot_all_sess_behav(Ax,hObject,type)

    handles = guidata(hObject) ;
    mouse = handles.mouse_menu.Value ;
    hh = handles.data_holder{mouse} ;
    
    Colors = {[0 0 1],[0.5 0.5 0.5],[1 0 0]} ;
    days = {'Baseline'} ;
    for i = 2 : numel(hh)
       days{i} = ['Bias ' num2str(i-1)] ;
    end
    
    for j = 1 : numel(hh)
        h = hh(j) ;
        trial_types = h.trial_types(h.ft:h.lt) ;
        trial_odors = h.trial_odors(h.ft:h.lt) ;

        for i = 1:3
            hit = sum(trial_types == 1 & trial_odors == i) ;
            CR = sum(trial_types == 2 & trial_odors == i) ;
            miss = sum(trial_types == 3 & trial_odors == i) ;
            FA = sum(trial_types == 4 & trial_odors == i) ;
            n_Go = hit+miss ;
            n_NoGo = CR+FA ;

            accuracy(i,j) = 100*(hit + CR)/(n_Go + n_NoGo) ;
            hit_rate(i,j) = 100*hit/n_Go ;
            FA_rate(i,j) = 100*FA/n_NoGo ;
            d_prime(i,j) = (norminv(min([hit_rate(i,j)/100 (n_Go-0.5)/n_Go])) - norminv(max([FA_rate(i,j)/100 0.5/n_NoGo]))) ;
            bias(i,j) = (norminv(min([hit_rate(i,j)/100 (n_Go-0.5)/n_Go])) + norminv(max([FA_rate(i,j)/100 0.5/n_NoGo]))) ;
            latency(i,j) = mean(h.trial_latencies(trial_types == 1 & trial_odors == i)) ;

        end
    end
    switch type
        case 'bias'
            behav = bias ;
            Title = 'Bias' ;
            YLIM = [-1 2] ;
        case 'latency'
            behav = latency ;
            Title = 'Latency (s)' ;
            YLIM = [-0.5 0.5] ;
        case 'FA rate'
            behav = FA_rate ;
            Title = 'FA rate' ;
            YLIM = [-50 50] ;
        case 'hit rate'
            behav = hit_rate ;
            Title = 'hit rate' ;
            YLIM = [-50 50] ;
    end
    plot(Ax,[0 (numel(hh)+1)],[0 0],'k')
    for i = [1 3]
        plot(Ax,behav(i,:)-behav(2,:),'Color',Colors{i})
    end
%     plot(Ax,bias(1,:)-bias(3,:),'k')    

    ylim(Ax , YLIM) 
    set(Ax,'XTick',1:numel(hh),'XTickLabel',days)
    ylabel(Ax,Title,'fontsize',16)
    
end

function plot_dprime(Ax,hObject)

    handles = guidata(hObject) ;
    legend(handles.anlss_ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    hit = [] ;
    miss = [] ;
    CR = [] ;
    FA = [] ;
    dprime = [] ;
    
    for i = h.ft:10:h.lt
        hit(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==1) ;
        CR(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==2) ;
        miss(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==3) ;
        FA(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==4) ;
        n = length(hit) ;
        m = max([1 n-9]) ;

        Go_lick_rate = sum(hit(m:n))/(sum(hit(m:n))+sum(miss(m:n))) ;
        NoGo_lick_rate = sum(FA(m:n))/(sum(FA(m:n))+sum(CR(m:n))) ;
                      
        dprime(end+1) = max(min(norminv(Go_lick_rate),2.3),-2.3) -...
            max(min(norminv(NoGo_lick_rate),2.3),-2.3) ;
    end 
    
    hold(Ax , 'on')
    plot(Ax , 1:length(dprime) , dprime , 'k') ;
    plot(Ax , 1:length(dprime) , dprime , '*k')
    ylim(Ax , [0 4.6]) 
    
    xlabel(Ax,'Trials/10','fontsize',16)
    ylabel(Ax,"d'",'fontsize',16)

end

function plot_lick_rate(Ax,hObject)

    handles = guidata(hObject) ;
    legend(handles.anlss_ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    trials_to_avarage = h.LICK_RATE_CRIT_TRIALS ;
    
    hit = [] ;
    miss = [] ;
    CR = [] ;
    FA = [] ;
    lick_rate = [] ;
    
    for i = h.ft:10:h.lt
        hit(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==1) ;
        CR(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==2) ;
        miss(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==3) ;
        FA(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==4) ;
        n = length(hit) ;
        m = max([n - trials_to_avarage/10 + 1 , 1]) ;
        lick_rate(end+1) = 100*(sum(hit(m:n)) + sum(FA(m:n)))/...
                          (sum(hit(m:n)) + sum(CR(m:n)) + sum(miss(m:n)) + sum(FA(m:n))) ;
    end 
    
    hold(Ax , 'on')
    plot(Ax , 1:length(lick_rate) , lick_rate , 'k') ;
    plot(Ax , 1:length(lick_rate) , lick_rate , '*k')
    ylim(Ax , [0 100]) 
    
    xlabel(Ax,'Trials/10','fontsize',16)
    ylabel(Ax,'lick rate (%)','fontsize',16)
    
end

function plot_raw_behav(hObject)

    n_trials = 30 ;

    handles = guidata(hObject) ;
    
    temp = figure() ;
    set(temp,'defaultAxesColorOrder',zeros(2,3))
    h_fig = gca ;
    hold on
    
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    p = h.Parameters ;
    trial_colors = {[0 0.7 0],[0 0.7 0],[0.7 0 0],[0.7 0 0]} ;
    texts1 = {'Hit' ,'CR' , 'Miss' , 'FA'} ;
    texts2 = {'Expected' ,'Neutral', 'Unexpected'} ;

    ft = h.ft ;
    types = h.trial_types ;
    odors = h.trial_odors ;
    tSOL = h.tSOL ;
    tSOS = h.tSOS ;
    tSOR = h.tSOR ;
    tone_dur = p.tone_dur*1e-3 ;
    odor_dur = p.odor_dur*1e-3 ;
    response_window = p.response_window ;
    
    count = n_trials ;
    i = 180 ;
    while count > 0
        if types(i) < 7
            rectangle(h_fig, 'Position' , [0 count-1 odor_dur 1], 'FaceColor' , handles.prob_colors{odors(i)})
            if types(i) < 5
                rectangle(h_fig,'Position' , [odor_dur count-1 tone_dur 1], 'FaceColor' , sum(types(i)==[1 3])*ones(1,3))
            else
                rectangle(h_fig,'Position' , [odor_dur count-1 tone_dur 1], 'FaceColor' , 'w')
            end
            rectangle(h_fig, 'Position' , [odor_dur+tone_dur count-1 response_window 1], 'FaceColor', 'w')
            licks = tSOL(tSOL < (tSOS(i)+tone_dur+response_window) & tSOL > tSOS(i)) - tSOS(i) + odor_dur ;
            scatter(h_fig , licks , (count-0.5)*ones(size(licks)) , 10 , 'k' ,'filled')
            rewards = tSOR(tSOR < (tSOS(i)+tone_dur+response_window+0.5) & tSOR > tSOS(i)) - tSOS(i) + odor_dur ;
            scatter(h_fig , rewards , (count-0.5)*ones(size(rewards)) , 40 , [0 0.69 0.235] ,'filled')
            if types(i) == 4
                scatter(h_fig , licks(1)+0.05 , count-0.5 , 40 , [0.753 0 0] ,'filled')
            end
            text(h_fig , 2.7 , count-0.5 , texts1{types(i)},'Color',trial_colors{types(i)})
            if (sum(types(i)==[1 3]) && odors(i) == 1) || (sum(types(i)==[2 4]) && odors(i) == 3)
                text(h_fig , 3 , count-0.5 , texts2{1},'Color',handles.exp_colors{1})
            elseif (sum(types(i)==[2 4]) && odors(i) == 1) || (sum(types(i)==[1 3]) && odors(i) == 3)
                text(h_fig , 3 , count-0.5 , texts2{3},'Color',handles.exp_colors{3})
            else
                text(h_fig , 3 , count-0.5 , texts2{2},'Color',handles.exp_colors{2})
            end
            count = count - 1 ;               
        end
        i = i + 1 ;
    end 
    
    axis(h_fig,[0 3.5 0 n_trials+0.5])
%     xlabel(h_fig,'Time (seconds)','fontsize',16)
%     ylabel(h_fig,'Trials','fontsize',16)
    set(h_fig,'YDir','reverse')
            
end

function plot_single(hObject,type)

    handles = guidata(hObject) ;
        
    figure()
    h_ax1 = gca ;
    hold on
    figure()
    h_ax2 = gca ;
    hold on

    if strcmp(handles.neuron_menu.String{handles.neuron_menu.Value},'No neurons')
        return
    end
    
    neur = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;    
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    n = h.Neurons(neur) ;
    
    odors = handles.odor_menu.Value - 1 ;
    if odors == 0
        odors = 1 : 3 ;
    elseif odors == 4
        odors = [1 3] ;        
    end

    switch handles.odor_types_menu.String{handles.odor_types_menu.Value}
        case 'Odor identities'
            baseline_colors = {[0 0.5 0] , [0.64 0.16 0.16] , [1 1 0]} ;
            colors = baseline_colors(h.odor_order(1:3)) ;
        case 'Odor bias'
            colors = handles.prob_colors ;
    end
    
    
    lims = [] ;
    n_frames = n.frames_before_stim + n.frames_for_amp ;
        

    for odor = odors
        dFF = n.get_dFF(odor,type) ;
        resps = n.get_responses(odor,type) ;
        if size(dFF,1) > 0
            x = h.Time ;
            switch handles.dFF_mode_menu.String{handles.dFF_mode_menu.Value}
                case 'Mean'
                    y = smoothdata(mean(dFF,1),'movmean',3) ;
                    dy = std(dFF,[],1)/sqrt(size(dFF,1)) ;
                    y = y(x>-1 & x<2.5) ;
                    dy = dy(x>-1 & x<2.5) ;
                    x = x(x>-1 & x<2.5) ;
                    plot(h_ax1,x, y , 'color' , colors{odor} , 'linewidth' , 4)
                    patch(h_ax1,[x flip(x)] , [y+dy flip(y-dy)]  , colors{odor} , 'facealpha' , 0.2)

                    lims(:,end+1) = [max(y(1:n_frames)) ; min(y(1:n_frames))] ;
                case 'Trials'
                    for trial = 1 : size(dFF,1)
                        y = dFF(trial,:) ;
                        plot(x, y , 'color' , colors{odor} , 'linewidth' , 1)
                        lims(:,end+1) = [max(y(1:n_frames)) ; min(y(1:n_frames))] ;
                    end
            end
            bar(h_ax2,find(odors==odor),mean(resps) ,'FaceColor', colors{odor})
            errorbar(h_ax2,find(odors==odor),mean(resps),std(resps)/length(resps),'.','Color','k','LineWidth',5,'CapSize',15)
        end
    end
    
    set(h_ax1,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')
    set(h_ax2,'XTick',[])
    plot(h_ax1,[-0.75 -0.75],[0.3 0.4]-0.2,'k','linewidth',2)
    plot(h_ax1,[-0.75 0.25],[0.3 0.3]-0.2,'k','linewidth',2)
    ht = text(h_ax1,-0.9,0.3-0.2,'0.1 dF/F','fontsize',14) ;
    text(h_ax1,-0.75,0.285-0.2,'1 second','fontsize',14)
    set(ht , 'rotation' , 90)
    rectangle(h_ax1,'Position',[-0.5 -0.05 0.5 0.01],'FaceColor','g')
    rectangle(h_ax1,'Position',[0 -0.05 0.15 0.01],'FaceColor',(type==1)*ones(1,3))


%     axis([-1 4 min(lims(2,:))-0.1 max(lims(1,:))+0.1]) ;
    axis(h_ax1,[-1 2.5 -0.15 0.255]) ;
    axis(h_ax2,[0 numel(odors)+1 0 0.3]) ;
    
end

function plot_single_sess_behav_all(Ax,hObject)

    handles = guidata(hObject) ;
    legend(Ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ; 
    
    hold(handles.anlss_ax , 'on')
    
    ylim2 = [-0.3 3] ;
    ylim1 = ylim2*110/ylim2(2) ;
    colororder(Ax,{'k','k'})
    
    hit_rate = h.hit_rate ;
    FA_rate = h.FA_rate ;
    bias = h.bias ;
    d_prime = h.d_prime_exp ;
    
    for odor = 1:3        
        bar(Ax, [0 4] + odor ,[hit_rate(odor) FA_rate(odor)],'FaceColor' , handles.prob_colors{odor},'BarWidth',0.1)
    end 
    
    ylim(Ax , ylim1)
    ylabel(Ax,'Rate (%)','fontsize',16)
    
    plot(Ax,[8 8],ylim1,'--k','Linewidth',2)
    
    yyaxis(Ax,'right')
        
    for odor = 1:3
        
        bar(Ax, 8 + odor ,bias(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.5)
        bar(Ax, 12 + odor ,d_prime(odor),'FaceColor' , handles.exp_colors{odor},'BarWidth',0.5)
    end 

    ylim(Ax , ylim2)
    ylabel(Ax,"d'/Bias (AU)",'fontsize',16)
     
    set(Ax,'XTick',(0:4:16) + 2,'XTickLabel',{'Hit rate','FA rate','Bias',"d'"})
    
end

function plot_single_sess_behav_parts(Ax,hObject,type)

    handles = guidata(hObject) ;
    legend(Ax,'hide')
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ; 
    
    hold(handles.anlss_ax , 'on')
        
    for odor = 1 : 3
        switch type
            case 'HR'
                hit_rate = h.hit_rate ;
                bar(Ax, odor ,hit_rate(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 100] ;
                Title = 'Hit rate (%)' ;
                yticks = [0 50 100] ;
                ylabels = {'0','50','100'} ;
            case 'FAR'
                FA_rate = h.FA_rate ;
                bar(Ax, odor ,FA_rate(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 100] ;
                Title = 'FA Rate (%)' ;
                yticks = [0 50 100] ;
                ylabels = {'0','50','100'} ;
            case 'bias'
                bias = h.bias ;
                bar(Ax, odor ,bias(odor),'FaceColor' , handles.prob_colors{odor},'BarWidth',0.8)
                YLIM = [0 2] ;
                Title = 'Bias' ;
                yticks = [0 1 2] ;
                ylabels = {'0','001','002'} ;
            case "d' expectation"
                d_prime = h.d_prime_exp ;
                bar(Ax, odor ,d_prime(odor),'FaceColor' , handles.exp_colors{odor},'BarWidth',0.8)
                YLIM = [0 2] ;
                Title = "d'" ;
                yticks = [0 1 2] ;
                ylabels = {'0','001','002'} ;
        end
    end
        
    
    axis(Ax , [0.3 3.7 YLIM])
%     ylabel(Ax,Title,'fontsize',16)     
    set(Ax,'XTick',[],'YTick',yticks,'YTickLabel',ylabels)
    
end

function plot_trial_counts(Ax,hObject)
    
    handles = guidata(hObject) ;
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    
    hit = [] ;
    miss = [] ;
    CR = [] ;
    FA = [] ;
    
    for i = h.ft:10:h.lt
        hit(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==1) ;
        CR(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==2) ;
        miss(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==3) ;
        FA(end+1) = sum(h.trial_types(i : min([i+10 length(h.trial_types)]))==4) ;        
    end
    
    hit = cumsum(hit) ;
    miss = cumsum(miss) ;
    CR = cumsum(CR) ;
    FA = cumsum(FA) ;
    
    hold(handles.anlss_ax , 'on')
    x(1) = plot(Ax , 1:length(hit) , hit , 'b' , 'DisplayName' , 'Hit') ;
    plot(Ax , 1:length(hit) , hit , '*b')
    x(2) = plot(Ax , 1:length(miss) , miss , 'r' , 'DisplayName' , 'Miss') ;
    plot(Ax , 1:length(miss) , miss , '*r')
    x(3) = plot(Ax , 1:length(CR) , CR , 'c' , 'DisplayName' , 'CR') ;
    plot(Ax , 1:length(CR) , CR , '*c')
    x(4) = plot(Ax , 1:length(FA) , FA , 'Color' , [0.8438 0.3281 0.0977] , 'DisplayName' , 'FA') ;
    plot(Ax , 1:length(FA) , FA , '*' , 'Color' , [0.8438 0.3281 0.0977])
    
    xlabel(Ax,'Trials/10','fontsize',16)
    ylabel(Ax,'#','fontsize',16)
    lgd = legend(x(1:4) , 'Location' , 'northwest') ; 
    lgd.FontSize = 16 ;
    
end

function plot_trial_sequence(hObject)

    handles = guidata(hObject) ;
    
    first_trial = 120 ;
    n_trials = 7 ;
    y_lims = [-0.1 0.5] ;
        
    figure()
    set(gcf,'position',[0 400 1700 300])
    h_ax = gca ;
    hold on

    if strcmp(handles.neuron_menu.String{handles.neuron_menu.Value},'No neurons')
        return
    end
    
    neur = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;    
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    n = h.Neurons(neur) ;

    odor_colors = handles.prob_colors ;    
    
    frames = h.frames_before_stim + (-10:(h.frames_for_amp+7)) ;
    x = h.Time(frames) ;
        
    types = h.trial_types ;
    odors = h.trial_odors ;
    for odor = 1 : 3
        for type = 1 : 4
            pos = types==type & odors==odor ;
            temp = n.valid_trials{odor,type} ;
            temp = type*temp ;
            temp(temp==0) = 7 ;
            types(pos) = temp ;
        end
    end
    
    plot(h_ax,[-3.5 -3.5],[0 0.2],'k','linewidth',2)
    plot(h_ax,[-3.5 -2.5],[0 0],'k','linewidth',2)
    ht = text(-3.8,0.022,'0.2 dF/F','fontsize',14) ;
    text(-3.5,-0.02,'1 sec','fontsize',14)
    set(ht , 'rotation' , 90)
    
    count = 0 ;
    ind = 0 ;
    x_start = 0 ;
    while count < n_trials
        trial = first_trial + ind ;
        ind = ind + 1 ;
        type = types(trial) ;
        if type~=7
            odor = odors(trial) ;
            dFF = n.get_dFF(odor,type) ;
            y = smoothdata(dFF(sum(types(1:trial)==type & odors(1:trial)==odor),:),'movmean',3) ;
            y = y(frames) ;
            plot(h_ax , x_start + x , y , 'k' , 'LineWidth' , 2)
            rectangle(h_ax,'Position',[x_start-0.5 -0.05 0.5 0.01],'FaceColor',odor_colors{odor})
            rectangle(h_ax,'Position',[x_start -0.05 0.15 0.01],'FaceColor',(mod(type,2)==1)*ones(1,3))
            x_start = x_start + (x(end)-x(1)) + 0.5 ;
            count = count + 1 ;
        end
    end 
    
    set(gca,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')

    ylim(y_lims) ;
    
end

function plot_all_neurons(hObject)

    handles = guidata(hObject) ;
        
    figure()
    h_ax = gca ;
    hold on
    
    Colors = handles.prob_colors ;

    mouse = handles.mouse_menu.Value ;
    sessions = [1 3] ;
    n1 = handles.data_holder{mouse}(sessions(1)).Neurons ;
    n2 = handles.data_holder{mouse}(sessions(2)).Neurons ;
    
    Time = handles.data_holder{mouse}(sessions(1)).Time ; 
    
    neurons = find(([n1.responsive] | [n2.responsive]) & (~[n1.invalid] & ~[n2.invalid])) ;
    neurons(22) = [] ;
%     neurons(end+1) = neurons(31) ;
%     neurons(31) = [] ;
    y = cell(length(neurons),12) ;
    dy = cell(length(neurons),12) ;
    mins = zeros(length(neurons),12) ;
    maxs = zeros(length(neurons),12) ;
    sig = zeros(length(neurons),2) ;
    for neur = 1 : length(neurons)
        for type = 1 : 2
            for odor = 1 : 3
                sig(neur,1) = n1(neurons(neur)).responsive ;
                sig(neur,2) = n2(neurons(neur)).responsive ;
                dFF1 = n1(neurons(neur)).get_dFF(odor,type) ;
                temp1 = smoothdata(mean(dFF1,1),'movmean',3) ;
                temp1 = temp1(Time>-1 & Time<2.5) ;
                y{neur,(type-1)*3+odor} = temp1 ;
                temp2 = std(dFF1,[],1)/sqrt(size(dFF1,1)) ;
                temp2 = temp2(Time>-1 & Time<2.5) ;
                dy{neur,(type-1)*3+odor} = temp2 ;
                mins(neur,(type-1)*3+odor) = min(temp1-temp2) ;
                maxs(neur,(type-1)*3+odor) = max(temp1+temp2) ;
                
                dFF2 = n2(neurons(neur)).get_dFF(odor,type) ;
                temp1 = smoothdata(mean(dFF2,1),'movmean',3) ;
                temp1 = temp1(Time>-1 & Time<2.5) ;
                y{neur,6+(type-1)*3+odor} = temp1 ;
                temp2 = std(dFF2,[],1)/sqrt(size(dFF2,1)) ;
                temp2 = temp2(Time>-1 & Time<2.5) ;
                dy{neur,6+(type-1)*3+odor} = temp2 ;
                mins(neur,6+(type-1)*3+odor) = min(temp1-temp2) ;
                maxs(neur,6+(type-1)*3+odor) = max(temp1+temp2) ;
            end
        end       
    end
    
    Time = Time(Time>-1 & Time<2.5) ;
    si = mode(diff(Time)) ;
    Max = max(maxs(:)) ;
    Min = min(mins(:)) ;
    Space = 0 ;
    for neur = 1 : size(y,1)
        for trace = 1 : 12
            yy = y{neur,trace} + sum(Space) ;
            dyy = dy{neur,trace} ;
            x = (0 : (length(yy)-1)) + (trace-1)*(length(yy)+1) + length(yy)*(trace>6) ;
            plot(h_ax , x , yy , 'color' , Colors{mod(trace-1,3)+1})
            patch(h_ax , [x flip(x)] ,[yy+dyy flip(yy-dyy)]  , Colors{mod(trace-1,3)+1} , 'facealpha' , 0.2)
%             patch(h_ax , [x flip(x)] ,[yy+dyy flip(yy-dyy)]  , Colors{mod(trace-1,3)+1})

%             plot(h_ax, ((neur-1)*(Max-Min) - 0.1)*[1 1] , 'LineWidth' , 4)
            

        end
        if sig(neur,1)
            scatter(6.4*length(Time),sum(Space),5,'k','filled')
        end
        if sig(neur,2)
            scatter(6.7*length(Time),sum(Space),5,'k','filled')
        end
        if neur < size(y,1)
            Space = [Space max(maxs(neur,:)) - min(mins(neur+1,:))] ;
        end
    end
    
    for line = 1 : 12
        plot(h_ax,find(Time==0)*[1 1] + (line-1)*(length(Time)+1) + length(Time)*(line>6), [-1 18],'color',[0.8 0.8 0.8])
    end
    
    plot(h_ax,[-3 -3],[-0.3 0.2],'k')
    plot(h_ax,[-3 -3+1/si],[-0.3 -0.3],'k')
    
    axis([-5 350 -0.5 11.5])
    
    set(h_ax,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')


    
end

% ========================== All mice analysis ===========================

function all_auditory_neurons(Ax,hObject,anlss_type,time,view)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    y_values = [] ;
    x_values = [] ;
    ids = [] ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    
    switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
        case 'Go neurons'
            pref_sound = 1 ;
        case 'NoGo neurons'
            pref_sound = 0 ;
        case 'None-selective'
            pref_sound = -1 ;
        otherwise
            pref_sound = [-1 0 1] ;
    end
    
    for mouse = mice
        switch time
            case 'this session'
                sessions = handles.session_menu.Value ;
            case 'bias'
                sessions = 3:5 ;
        end
        for sess = sessions
            hh = h{mouse}(sess) ;
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive && sum(n.pref_sound == pref_sound)
                    switch anlss_type
                        case 'expectation DI'
                            y_values = [y_values ; n.exp_DI] ;
                            xloc = 1:3 ;
                            x_values = [x_values ; xloc] ;
                            YLIM = [0.4 1.2] ;
                            yticks = 0.5:0.1:1 ;
                            Title = 'DI' ;
                            XLabels = {'Expected','Neutral','Unexpected'} ;
                        case 'expectation DI diff'
                            y_values = [y_values ; n.exp_DI_diff] ;
                            xloc = 1 ;
                            x_values = [x_values ; xloc] ;
                            YLIM = [-0.55 0.55] ;
                            yticks = -0.5:0.1:0.5 ;
                            Title = '\DeltaDI' ;
                            XLabels = {'Unexpected-Expected'} ;
                    end
                    ids(end+1,:) = [mouse sess neur] ;
                end
            end
        end
    end
    
    handles.h_scat = {} ;
    for neur = 1 : size(y_values,1)
        handles.h_scat{neur} = scatter(Ax,...
            x_values(neur,:)-0.2+0.4*rand(size(x_values(neur,:))),y_values(neur,:),20,...
            'MarkerEdgeColor','k','MarkerFaceColor','k') ;
        set(handles.h_scat{neur},'ButtonDownFcn',{@scatter_func , hObject,ids(neur,:),2}) ;
    end

    for i = 1 : 3
        Data{i} = y_values(x_values==i) ;
        plot(Ax,i,mean(y_values(x_values==i)),'+r','markersize',20)
    end
    
    if ~strcmp(anlss_type,'expectation DI diff')
        p(1) = signrank(Data{1},Data{2}) ;
        text(Ax,1.2,1.05,sig_signs{sum(3*p(1)<=pp)},'fontsize',16)
        plot(Ax,[1. 1.9],[1 1],'k')
        p(2) = signrank(Data{2},Data{3}) ;
        text(Ax,2.2,1.05,sig_signs{sum(3*p(2)<=pp)},'fontsize',16)
        plot(Ax,[2.1 2.9],[1 1],'k')
        p(3) = signrank(Data{1},Data{3}) ;
        text(Ax,1.7,1.15,sig_signs{sum(3*p(3)<=pp)},'fontsize',16)
        plot(Ax,[1.1 2.9],[1.1 1.1],'k')
    end

    axis(Ax , [0 4 YLIM]) 
    set(Ax,'XTick',xloc,'XTickLabel',XLabels,'YTick',yticks,'fontsize',16)
    ylabel(Ax,Title,'fontsize',16)
    
    if strcmp(view,'violin')
        figure
        boxplot(y_values)
        hold on
        text(gca,1.2,1.05,sig_signs{sum(3*p(1)<=pp)},'fontsize',16)
        plot(gca,[1. 1.9],[1 1],'k')
        text(gca,2.2,1.05,sig_signs{sum(3*p(2)<=pp)},'fontsize',16)
        plot(gca,[2.1 2.9],[1 1],'k')
        text(gca,1.7,1.15,sig_signs{sum(3*p(3)<=pp)},'fontsize',16)
        plot(gca,[1.1 2.9],[1.1 1.1],'k')
        axis(gca , [0 4 YLIM]) 
        set(gca,'XTick',[],'YTick',yticks,'fontsize',16)
%         ylabel(gca,Title,'fontsize',16)
    end
    
    guidata(hObject,handles)

end

function plot_odor_behav_prog(Ax,hObject,type)

    handles = guidata(hObject) ;
    
    Colors = handles.prob_colors ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [1,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
    
    behav = cell(1,max(n_days)) ;
    
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : n_days(mouse)
            h = hh(day) ;               
            switch type
                case 'bias'
                    behav{day}(:,end+1) = h.bias ;
                    Title = '\DeltaBias' ;
                    YLIM = [-1 1] ;
                    yticks = [YLIM(1) 0 YLIM(2)] ;
                    ylabels = {'-01','0','01'} ;
                case 'latency'
                    behav{day}(:,end+1) = h.hit_latency ;
                    Title = '\DeltaLatency (sec)' ;
                    YLIM = [-0.1 0.4] ;
                case 'FA rate'
                    behav{day}(:,end+1) = h.FA_rate(1:3) ;
                    Title = '\DeltaFA rate' ;
                    YLIM = [-30 30] ;
                    yticks = [YLIM(1) 0 YLIM(2)] ;
                    ylabels = {num2str(YLIM(1)),'0',num2str(YLIM(2))} ;
                case 'hit rate'
                    behav{day}(:,end+1) = h.hit_rate(1:3) ;
                    Title = '\DeltaHit rate' ;
                    YLIM = [-30 30] ;
                    yticks = [YLIM(1) 0 YLIM(2)] ;
                    ylabels = {num2str(YLIM(1)),'0',num2str(YLIM(2))} ;
                case "d' expectation"
                    behav{day}(:,end+1) = h.d_prime_exp ;
                    Title = "\Deltad'" ;
                    YLIM = [-1 1] ;
                    yticks = [YLIM(1) 0 YLIM(2)] ;
                    ylabels = {'-01','0','01'} ;
                    Colors = handles.exp_colors ;

            end
        end
    end
    
    for day = 1 : numel(behav)
        d_behav{day} = behav{day} - repmat(behav{day}(2,:),3,1) ;
        m_behav(:,day) = mean(d_behav{day},2) ;
        s_behav(:,day) = std(d_behav{day},[],2)/sqrt(size(behav{day},2)) ;
    end
    
    for odor = [1 3]
        errorbar(Ax,m_behav(odor,:),s_behav(odor,:),'Color',Colors{odor},'LineWidth',2)
    end
    
    plot(Ax,0:numel(behav)+1,zeros(1,numel(behav)+2),'k')
    
    for day = 1 : numel(d_behav)
        odor1 = d_behav{day}(1,:) ;
        odor2 = d_behav{day}(2,:) ;
        odor3 = d_behav{day}(3,:) ;
        
        p = signrank(odor1,odor3,'Tail','Right') ;
        if ~ isnan(p)
            text(Ax,day,0.9*YLIM(2),sig_signs{sum(p<=pp)},'fontsize',16)
        end
        
    end
    

    days = {'PBS'} ;
    for i = 1 : numel(behav)
        days{i+1} = ['BS ' num2str(i)] ;
    end
    axis(Ax , [0.5 7.5 YLIM]) 
    set(Ax,'XTick',1:numel(behav),'XTickLabel',{},'YTick',yticks,'YTickLabel',ylabels)
%     ylabel(Ax,Title,'fontsize',16)
    
end

function plot_odor_sound_behav_prog(Ax,hObject)

    handles = guidata(hObject) ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [1,0.05,0.005,0.0005,0.00005,0.000005] ;
    YLIM = [-1 2.5] ;
       
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
       
    odor_dprime = cell(1,max(n_days)) ;
    sound_dprime = cell(1,max(n_days)) ;    
        
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : numel(hh)           
            h = hh(day) ;

            sound_dprime{day}(end+1) = h.d_prime(4) ;            
            odor_dprime{day}(end+1) = h.odor_d_prime ;
        end
    end
    
 
    for day = 1 : numel(sound_dprime)
        m_sound_behav(day) = mean(sound_dprime{day}) ;
        s_sound_behav(day) = std(sound_dprime{day})/sqrt(numel(sound_dprime{day})) ;
        m_odor_behav(day) = mean(odor_dprime{day}) ;
        s_odor_behav(day) = std(odor_dprime{day})/sqrt(numel(odor_dprime{day})) ;
    end
   
    errorbar(Ax,m_sound_behav,s_sound_behav,'k') 
    errorbar(Ax,m_odor_behav,s_odor_behav,'--k')
    plot(Ax,[0 numel(m_sound_behav)+1],[0 0],'k')
    
    m_sound_behav
    s_sound_behav
    m_odor_behav
    s_odor_behav
    
    for day = 1 : numel(m_sound_behav)
        
        p(day) = signrank(sound_dprime{day},odor_dprime{day},'Tail','Right') ;
        if ~ isnan(p(day))
            text(Ax,day,0.9*YLIM(2),sig_signs{sum(p(day)<=pp)},'fontsize',16)
        end
        
    end
    

%     days = {'PBS'} ;
%     for i = 1 : numel(m_sound_behav)
%         days{i+1} = ['BS ' num2str(i)] ;
%     end
%     axis(Ax , [0 max(n_days)+1 YLIM]) 
%     set(Ax,'XTick',1:numel(m_sound_behav),'XTickLabel',days)
    set(Ax,'XTick',1:numel(m_sound_behav),'XTickLabel',{})
%     ylabel(Ax,"d'",'fontsize',16)
    
end

function plot_select_days_behav(Ax,hObject,type)

    handles = guidata(hObject) ;
    
    first_day = 3 ;
    last_day = 5 ;
    
    
    Colors = handles.prob_colors ;
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    for mouse = 1 : numel(handles.data_holder)
        n_days(mouse) = numel(handles.data_holder{mouse}) ;
    end
    
    behav = cell(1,max(n_days)) ;
    
    for mouse = 1 : numel(handles.data_holder)
        hh = handles.data_holder{mouse} ;
        for day = 1 : n_days(mouse)
            h = hh(day) ;

            switch type
                case 'bias'
                    behav{day}(:,end+1) = h.bias(1:3) ;
                    Title = 'Bias' ;
                    YLIM = [-2 5] ;
                    yticks = -2:2:4 ;
                    ylabels = {'-02','000','002','004'} ;
                case 'latency'
                    behav{day}(:,end+1) = h.hit_latency(1:3) ;
                    Title = 'Latency (sec)' ;
                    YLIM = [-0.1 0.4] ;
                case 'FAR'
                    behav{day}(:,end+1) = h.FA_rate(1:3) ;
                    Title = 'FA rate' ;
                    YLIM = [0 110] ;
                    yticks = [0 50 100] ;
                    ylabels = {'0','50','100'} ;
                case 'HR'
                    behav{day}(:,end+1) = h.hit_rate(1:3) ;
                    Title = 'Hit rate' ;
                    YLIM = [0 110] ;
                    yticks = [0 50 100] ;
                    ylabels = {'0','50','100'} ;
                case "d' expectation"
                    behav{day}(:,end+1) = h.d_prime_exp(1:3) ;
                    Title = "d'" ;
                    YLIM = [0 4] ;
                    yticks = 0:2:4 ;
                    ylabels = {'000','002','004'} ;
                    Colors = handles.exp_colors ;

            end
        end
    end
       
    d_behav = [] ;
    for day = first_day : last_day
        d_behav = [d_behav behav{day}] ;
    end
    m_behav = mean(d_behav,2) ;
    s_behav = std(d_behav,[],2)/sqrt(size(behav,2)) ;
    
    for odor = 1:3
        scatter(Ax,odor*ones(1,size(d_behav,2)),d_behav(odor,:),100,Colors{odor},'filled')
%         errorbar(Ax,odor,m_behav(odor),s_behav(odor),'.','Color','k','LineWidth',7,'CapSize',15)
    end
    for sess = 1 : size(d_behav,2)
        plot(1:3,d_behav(:,sess),'k')
    end
    
%     [~,p] = ttest(d_behav(1,:),d_behav(3,:),'Tail','Right') ;
%     if ~ isnan(p)
%         text(Ax,2,0.9*YLIM(2),sig_signs{sum(3*p<=pp)},'fontsize',16)
%         plot(Ax,[1 3],0.8*YLIM(2)*[1 1],'k')
%     end
% 
%     [~,p] = ttest(d_behav(1,:),d_behav(2,:),'Tail','Right') ;
%     if ~ isnan(p)
%         text(Ax,1.5,0.7*YLIM(2), sig_signs{sum(3*p<=pp)},'fontsize',16)
%         plot(Ax,[1.1 1.9],0.6*YLIM(2)*[1 1],'k')
%     end
% 
%     [~,p] = ttest(d_behav(3,:),d_behav(2,:),'Tail','Left') ;
%     if ~ isnan(p)
%         text(Ax,2.5,0.7*YLIM(2), sig_signs{sum(3*p<=pp)},'fontsize',16)
%         plot(Ax,[2.1 2.9],0.6*YLIM(2)*[1 1],'k')
%     end
    


    axis(Ax , [0.5 3.5 YLIM]) 
    set(Ax,'XTick',[],'YTick',yticks,'YTickLabel',ylabels)
%     ylabel(Ax,Title,'fontsize',16)
    
end

function plot_DI_progress(Ax,hObject)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    sessions = 1 : 7 ;
    DI = cell(1,length(sessions)) ;
    
    switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
        case 'Go neurons'
            pref_sound = 1 ;
        case 'NoGo neurons'
            pref_sound = 0 ;
        case 'None-selective'
            pref_sound = -1 ;
        otherwise
            pref_sound = [-1 0 1] ;
    end
    
    for mouse = mice
        for sess = sessions
            hh = h{mouse}(sess) ; 
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive && sum(n.pref_sound == pref_sound)
                    DI{sess} = [DI{sess} ; n.exp_DI([1 3])] ;
                end
            end
        end
    end
    
    for sess = sessions
        m_DI(sess,:) = mean(DI{sess},1) ;
        s_DI(sess,:) = std(DI{sess},[],1)/sqrt(size(DI{sess},1)) ;
        p = signrank(DI{sess}(:,1),DI{sess}(:,2)) ;
        p_text{sess} = sig_signs{sum(p<=pp)} ;
    end
    errorbar(Ax,m_DI(:,1),s_DI(:,1),'Color',handles.exp_colors{1},'LineWidth',2)
    errorbar(Ax,m_DI(:,2),s_DI(:,2),'Color',handles.exp_colors{3},'LineWidth',2)
    text(Ax,sessions,0.67*ones(size(sessions)),p_text,'fontsize',16)
    
    axis(Ax,[0 length(sessions)+1 0.55 0.7])
    days = {'PBS'} ;
    for i = 1 : numel(sessions)
        days{i+1} = ['BS ' num2str(i)] ;
    end
%     set(Ax,'Xtick',1:length(sessions),'XTickLabel',days)
    set(Ax,'Xtick',1:length(sessions),'XTickLabel',{},'YTick',[0.55 0.6 0.65 0.7])
%     ylabel(Ax,'DI','Fontsize',16)
        
        
    guidata(hObject,handles)
    
end

function plot_DI_diff_progress(Ax,hObject)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    
    mice = handles.phys_mice ;
    sessions = 1 : 7 ;
    dDI = cell(1,length(sessions)) ;
    
    switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
        case 'Go neurons'
            pref_sound = 1 ;
        case 'NoGo neurons'
            pref_sound = 0 ;
        case 'None-selective'
            pref_sound = -1 ;
        otherwise
            pref_sound = [-1 0 1] ;
    end
    
    for mouse = mice
        for sess = sessions
            hh = h{mouse}(sess) ; 
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive && sum(n.pref_sound == pref_sound)
                    dDI{sess}(end+1) = -diff(n.exp_DI([1 3])) ;
                end
            end
        end
    end
    
    for sess = sessions
        m_dDI(sess,:) = mean(dDI{sess}) ;
        s_dDI(sess,:) = std(dDI{sess})/sqrt(length(dDI{sess})) ;
        p = signrank(dDI{sess},0) ;
        p_text{sess} = sig_signs{sum(p<=pp)} ;
    end
    plot(Ax,[0 8],[0 0],'--k')
    errorbar(Ax,m_dDI,s_dDI,'Color','k','LineWidth',2)
    text(Ax,sessions,0.04*ones(size(sessions)),p_text,'fontsize',16)
    
    axis(Ax,[0 length(sessions)+1 -0.05 0.02])
    days = {'PBS'} ;
    for i = 1 : numel(sessions)
        days{i+1} = ['BS ' num2str(i)] ;
    end
%     set(Ax,'Xtick',1:length(sessions),'XTickLabel',days)
    set(Ax,'Xtick',1:length(sessions),'XTickLabel',{},'YTick',[-0.05 -0.025 0])
%     ylabel(Ax,'DI','Fontsize',16)
        
        
    guidata(hObject,handles)
    
end

function plot_all_population_means(hObject,time)

    plot_population_mean(hObject,time,1)
    plot_population_mean(hObject,time,0)

end

function plot_population_mean(hObject,time,neuron_type)

    handles = guidata(hObject) ;
    h = handles.data_holder ;
    
    if length(neuron_type)==1
        if neuron_type
            Title = 'Go neurons' ;
        else
            Title = 'NoGo neurons' ;
        end
    else 
        Title = 'All neurons' ;
    end
        
    sig_signs = {'n.s.','*','**','***','****','*****'} ;
    pp = [inf,0.05,0.005,0.0005,0.00005,0.000005] ;
    colors = handles.prob_colors ;
        
    mice = handles.phys_mice ;
    
    traces = cell(2,3) ;
    resps = [] ;
    
    temp = {} ;
    
    for mouse = mice
        switch time
            case 'this session'
                sessions = handles.session_menu.Value ;
            case 'bias'
                sessions = 3:5 ;
        end
        for sess = sessions
            hh = h{mouse}(sess) ;
            for neur = 1 : hh.n_neurons
                n = hh.Neurons(neur) ;
                if ~n.invalid && n.responsive && sum(n.pref_sound==neuron_type)
                    temp{end+1} = cell(2,3) ;
                    resps(end+1,:,:) = [n.mean_resps(:,1:2) mean(n.mean_resps(:,1:2),2)] ;
                    for type = 1:2
                        for odor = 1 : 3
                            traces{type,odor} = [traces{type,odor} ; mean(n.get_dFF(odor,type),1)] ;
                        end
                    end
                end
            end
        end
    end 
    
    figure()
    
    y_lim = 0 ;
    for odor = 1:3
        x = hh.Time ;
        y1 = mean(traces{1,odor},1) ;
        dy1 = std(traces{1,odor},[],1)/sqrt(size(traces{1,odor},1)) ;
        y2 = mean(traces{2,odor},1) ;
        dy2 = std(traces{2,odor},[],1)/sqrt(size(traces{2,odor},1)) ;
        y1 = y1(x>-1 & x<2) ;
        dy1 = dy1(x>-1 & x<2) ;
        y2 = y2(x>-1 & x<2) ;
        dy2 = dy2(x>-1 & x<2) ;
        x = x(x>-1 & x<2) ;
        
%         y_lim = max([y_lim 1.2*max(y1) 1.2*max(y2)]) ;
        
        subplot(3,2,[1 3])
        plot(x, y1 , 'color' , colors{odor} , 'linewidth' , 2)
        patch([x flip(x)] , [y1+dy1 flip(y1-dy1)]  , colors{odor} , 'facealpha' , 0.2)
        axis([-1 2 -0.04 0.11])
        hold on
        
        set(gca,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')
        plot(gca,[-0.6 -0.6],[0.07 0.09],'k','linewidth',2)
        plot(gca,[-0.6 0.4],[0.07 0.07],'k','linewidth',2)
        ht = text(gca,-0.85,0.07,'0.02 dF/F','fontsize',10) ;
        text(gca,-0.75,0.065,'1 sec','fontsize',10)
        set(ht , 'rotation' , 90)
        rectangle(gca,'Position',[-0.5 -0.015 0.5 0.005],'FaceColor','g')
        rectangle(gca,'Position',[0 -0.015 0.15 0.005],'FaceColor',(type==1)*ones(1,3))
        
        subplot(3,2,5)
        hold on 
        bar(odor,mean(resps(:,odor,1)),'FaceColor',colors{odor})
        errorbar(odor,mean(resps(:,odor,1)),std(resps(:,odor,1))/sqrt(size(resps,1)),'k')
        set(gca,'XTick',[])
        ylim([0 0.12])     
        
        subplot(3,2,[2 4])
%         title(Title,'FontSize',16)
        plot(x, y2 , 'color' , colors{odor} , 'linewidth' , 2)
        patch([x flip(x)] , [y2+dy2 flip(y2-dy2)]  , colors{odor} , 'facealpha' , 0.2)
        axis([-1 2 -0.04 0.11])
        hold on
        
        set(gca,'XTick',[],'YTick',[],'box', 'off','xcolor','none','ycolor','none')
        rectangle(gca,'Position',[-0.5 -0.015 0.5 0.005],'FaceColor','g')
        rectangle(gca,'Position',[0 -0.015 0.15 0.005],'FaceColor',(type==1)*ones(1,3))
        
        subplot(3,2,6)
        hold on
        bar(odor,mean(resps(:,odor,2)),'FaceColor',colors{odor})
        errorbar(odor,mean(resps(:,odor,2)),std(resps(:,odor,2))/sqrt(size(resps,1)),'k')
        set(gca,'XTick',[])
        ylim([0 0.12])
        
    end
    
    subplot(3,2,5)
    [~,p] = ttest(resps(:,1,1),resps(:,2,1)) ;
    text(1.3,0.088,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([0.9 1.9],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,2,1),resps(:,3,1)) ;
    text(2.3,0.088,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([2.1 3.1],0.077*[1 1],'k')
        
    [~,p] = ttest(resps(:,1,1),resps(:,3,1)) ;
    text(1.8,0.11,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([1 3],0.1*[1 1],'k')
    
    subplot(3,2,6)
    [~,p] = ttest(resps(:,1,2),resps(:,2,2)) ;
    text(1.3,0.088,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([0.9 1.9],0.077*[1 1],'k')
    
    [~,p] = ttest(resps(:,2,2),resps(:,3,2)) ;
    text(2.3,0.088,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([2.1 3.1],0.077*[1 1],'k')
        
    [~,p] = ttest(resps(:,1,2),resps(:,3,2)) ;
    text(1.8,0.11,sig_signs{sum(3*p<=pp)},'fontsize',12)
    plot([1 3],0.1*[1 1],'k')
    
%     figure()
%     hold on
%     scatter(resps(:,1,2),resps(:,3,2))
%     plot([-0.3 0.3],[-0.3 0.3])
%     sum(resps(:,3,2)>resps(:,1,2))/size(resps,1)
    
    
end

function plot_dDI_to_dprime(Ax,hObject)

    handles = guidata(hObject) ;
    hh = handles.data_holder ;
    
    lims = [-1 1 -0.1 0.1] ;
    
    mice = handles.phys_mice ;
    
    switch handles.neuron_filter_menu.String{handles.neuron_filter_menu.Value}
        case 'Go neurons'
            pref_sound = 2 ;
        case 'NoGo neurons'
            pref_sound = 3 ;
        otherwise
            pref_sound = 1 ;
    end
    
    d_d_prime = [] ;
    dDI = [] ;
    handles.h_scat = {} ;
    
    for mouse = 1 : length(mice)
        for sess = 1 : (numel(hh{mice(mouse)}))
            h = hh{mice(mouse)}(sess) ;
            if ~isnan(h.mean_dDI(1))
                d_d_prime(end+1) = (h.d_prime_exp(1) - h.d_prime_exp(3))/(h.d_prime_exp(1) + h.d_prime_exp(3)) ;
                dDI(end+1) = -h.mean_dDI(pref_sound) ;
                handles.h_scat{mouse}{sess} = scatter(Ax,d_d_prime(end),dDI(end),40,...
                    'MarkerEdgeColor','k','MarkerFaceColor','k') ;
                set(handles.h_scat{mouse}{sess},'ButtonDownFcn',{@scatter_func , hObject,[mice(mouse) sess 1],1})
            end
        end
    end
    
    
    [R,p] = corrcoef(d_d_prime,dDI)
    
    plot(Ax,lims(1:2),[0 0],'k')
    plot(Ax,[0 0],lims(3:4),'k')
    
    par = polyfit(d_d_prime,dDI,1)
    
    plot(Ax,[-1 1],polyval(par,[-1.45 1.45]) , 'r','LineWidth',2)
    
%     xlabel(Ax,"d' CI",'FontSize',16)
%     ylabel(Ax,"<\DeltaDI>",'FontSize',16)

    set(Ax,'XTick',[lims(1) 0 lims(2)],'YTick',[lims(3) 0 lims(4)])
    
    guidata(hObject,handles)

end

% ========================== Helper functions ============================

function plot_dFF(hObject)

    handles = guidata(hObject) ;
    
    delete(get(handles.neuron_panel,'children'))
    
    h_fig = handles.neuron_panel ;

    if strcmp(handles.neuron_menu.String{handles.neuron_menu.Value},'No neurons')
        return
    end
    
    neur = str2double(handles.neuron_menu.String{handles.neuron_menu.Value}(8:end)) ;    
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse}(sess) ;
    if isnan(neur)
        return
    else
        n = h.Neurons(neur) ;
    end
    
    odors = handles.odor_menu.Value - 1 ;
    if odors == 0
        odors = 1 : 3 ;
    elseif odors == 4
        odors = [1 3] ;        
    end

    titles = {'Hit' , 'CR' , 'Miss' , 'FA'} ;
    colors = handles.prob_colors ;
    
    
    lims = [] ;
    n_frames = h.frames_before_stim + h.frames_for_amp ;
        
    for type = 1 : 4
        subplot(2,2,type,'Parent',h_fig)
        hold on
        for odor = odors
            dFF = n.get_dFF(odor,type) ;
            if size(dFF,1) > 0
                x = h.Time ;
                switch handles.dFF_mode_menu.String{handles.dFF_mode_menu.Value}
                    case 'Mean'
                        y = smoothdata(mean(dFF,1),'movmean',3) ;

                        dy = std(dFF,[],1)/sqrt(size(dFF,1)) ;
                        plot(x, y , 'color' , colors{odor} , 'linewidth' , 2)
                        patch([x flip(x)] , [y+dy flip(y-dy)]  , colors{odor} , 'facealpha' , 0.2)

                        if type < 3
                            lims(:,end+1) = [max(y(1:n_frames)) ; min(y(1:n_frames))] ;
                        end
                    case 'Trials'
                        for trial = 1 : size(dFF,1)
                            y = dFF(trial,:) ;
                            plot(x, y , 'color' , colors{odor} , 'linewidth' , 1)
                            if type < 3
                                lims(:,end+1) = [max(y(1:n_frames)) ; min(y(1:n_frames))] ;
                            end
                        end
                end

            end
        end
        title(titles{type})
        if type == 1
            xlabel('Time (seconds)')
            ylabel('\DeltaF/F_0')
        end

    end

    for type = 1 : 2
        subplot(2,2,type,'Parent',h_fig)
        axis([-1 4 min(lims(2,:))-0.1 max(lims(1,:))+0.1]) ;
        if sum(n.responsiveness(:,type)>0)
            plot(-0.4 , max(lims(1,:)) , '*k')
        end
        if n.inhibited
            plot(-0.2 , max(lims(1,:)) , '*r')
        end
    end
        
    
end

function scatter_func(src,~,hObject,ids,type)

    handles = guidata(hObject) ;

    handles.mouse_menu.Value = ids(1) ;
    mouse_menu_Callback(hObject, 0 , handles)
    handles.session_menu.Value = ids(2) ;
    session_menu_Callback(hObject, 0 , handles)
    handles.neuron_filter_menu.Value = 1 ;
    neuron_filter_menu_Callback(hObject, 0 , handles)
    handles.neuron_menu.Value = ids(3) ;
    neuron_menu_Callback(hObject, 0 , handles)
    
    for i = 1 : numel(handles.h_scat)
        for j = 1 : numel(handles.h_scat{i})
            switch type
                case 1
                    set(handles.h_scat{i}{j},'MarkerFaceColor','k')
                case 2
                    set(handles.h_scat{i},'MarkerFaceColor','k')
            end
        end
    end
        
    src.MarkerFaceColor = 'g' ;    
    
end

function scroll_func(hObject,eventdata)

    handles = guidata(hObject) ;
    
    fig_pos = get(hObject,'Position') ;
    mouse_pos = get(hObject,'CurrentPoint') ;
    
    anlss_type_menu_pos(1) = fig_pos(3) * handles.all_mice_anlss_menu.Position(1) ;
    anlss_type_menu_pos(2) = fig_pos(3) * (handles.all_mice_anlss_menu.Position(1) + handles.all_mice_anlss_menu.Position(3)) ;
    anlss_type_menu_pos(3) = fig_pos(4) * handles.all_mice_anlss_menu.Position(2) ;
    anlss_type_menu_pos(4) = fig_pos(4) * (handles.all_mice_anlss_menu.Position(2) + handles.all_mice_anlss_menu.Position(4)) ;
    
    anlss_type_menu_pos2(1) = fig_pos(3) * handles.one_mouse_anlss_menu.Position(1) ;
    anlss_type_menu_pos2(2) = fig_pos(3) * (handles.one_mouse_anlss_menu.Position(1) + handles.one_mouse_anlss_menu.Position(3)) ;
    anlss_type_menu_pos2(3) = fig_pos(4) * handles.one_mouse_anlss_menu.Position(2) ;
    anlss_type_menu_pos2(4) = fig_pos(4) * (handles.one_mouse_anlss_menu.Position(2) + handles.one_mouse_anlss_menu.Position(4)) ;
    
    session_menu_pos(1) = fig_pos(3) * handles.session_menu.Position(1) ;
    session_menu_pos(2) = fig_pos(3) * (handles.session_menu.Position(1) + handles.session_menu.Position(3)) ;
    session_menu_pos(3) = fig_pos(4) * handles.session_menu.Position(2) ;
    session_menu_pos(4) = fig_pos(4) * (handles.session_menu.Position(2) + handles.session_menu.Position(4)) ;
    
    neuron_menu_pos(1) = fig_pos(3) * handles.neuron_menu.Position(1) ;
    neuron_menu_pos(2) = fig_pos(3) * (handles.neuron_menu.Position(1) + handles.neuron_menu.Position(3)) ;
    neuron_menu_pos(3) = fig_pos(4) * handles.neuron_menu.Position(2) ;
    neuron_menu_pos(4) = fig_pos(4) * (handles.neuron_menu.Position(2) + handles.neuron_menu.Position(4)) ;
    
    neuron_filter_menu_pos(1) = fig_pos(3) * handles.neuron_filter_menu.Position(1) ;
    neuron_filter_menu_pos(2) = fig_pos(3) * (handles.neuron_filter_menu.Position(1) + handles.neuron_filter_menu.Position(3)) ;
    neuron_filter_menu_pos(3) = fig_pos(4) * handles.neuron_filter_menu.Position(2) ;
    neuron_filter_menu_pos(4) = fig_pos(4) * (handles.neuron_filter_menu.Position(2) + handles.neuron_filter_menu.Position(4)) ;
    
    if (mouse_pos(1) > anlss_type_menu_pos(1)) && (mouse_pos(1) < anlss_type_menu_pos(2)) &&...
       (mouse_pos(2) > anlss_type_menu_pos(3)) && (mouse_pos(2) < anlss_type_menu_pos(4))
   
        val = handles.all_mice_anlss_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.all_mice_anlss_menu.String) ;
   
        handles.all_mice_anlss_menu.Value = min([max([new_val 1]) max_val]) ;
        
        all_mice_anlss_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > anlss_type_menu_pos2(1)) && (mouse_pos(1) < anlss_type_menu_pos2(2)) &&...
       (mouse_pos(2) > anlss_type_menu_pos2(3)) && (mouse_pos(2) < anlss_type_menu_pos2(4))
   
        val = handles.one_mouse_anlss_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.one_mouse_anlss_menu.String) ;
   
        handles.one_mouse_anlss_menu.Value = min([max([new_val 1]) max_val]) ;
        
        one_mouse_anlss_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > session_menu_pos(1)) && (mouse_pos(1) < session_menu_pos(2)) &&...
           (mouse_pos(2) > session_menu_pos(3)) && (mouse_pos(2) < session_menu_pos(4))
   
        val = handles.session_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.session_menu.String) ;
   
        handles.session_menu.Value = min([max([new_val 1]) max_val]) ;
        
        session_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > neuron_menu_pos(1)) && (mouse_pos(1) < neuron_menu_pos(2)) &&...
           (mouse_pos(2) > neuron_menu_pos(3)) && (mouse_pos(2) < neuron_menu_pos(4))
   
        val = handles.neuron_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.neuron_menu.String) ;
   
        handles.neuron_menu.Value = min([max([new_val 1]) max_val]) ;
        
        neuron_menu_Callback(hObject,0,handles)
        
    elseif (mouse_pos(1) > neuron_filter_menu_pos(1)) && (mouse_pos(1) < neuron_filter_menu_pos(2)) &&...
           (mouse_pos(2) > neuron_filter_menu_pos(3)) && (mouse_pos(2) < neuron_filter_menu_pos(4))
   
        val = handles.neuron_filter_menu.Value ;
        new_val = val + eventdata.VerticalScrollCount ;
        max_val = numel(handles.neuron_filter_menu.String) ;
   
        handles.neuron_filter_menu.Value = min([max([new_val 1]) max_val]) ;
        
        neuron_filter_menu_Callback(hObject,0,handles)
               
    end
    
end

function update_probs_panel(hObject)

    handles = guidata(hObject) ;
    delete(get(handles.probs_panel,'children'))
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    p = handles.data_holder{mouse}(sess).Parameters ;
    
    ax = axes(handles.probs_panel) ;
    
    probs = [p.Go_odor1_prob/(p.Go_odor1_prob+p.NoGo_odor1_prob),...
             p.Go_odor2_prob/(p.Go_odor2_prob+p.NoGo_odor2_prob),...
             p.Go_odor3_prob/(p.Go_odor3_prob+p.NoGo_odor3_prob)] ;
    
    bar(ax,1,probs(1),'FaceColor','k')
    hold on
    bar(ax,2,probs(2),'FaceColor','k')
    bar(ax,3,probs(3),'FaceColor','k')
    set(ax,'XTick',1:3,'XTickLabel',{'EB','\alphaP','IAA'},'FontSize',10)
    ylim(ax,[0 1.1])
    ylabel('Go prob.','FontSize',10)    
    
end

function update_session_menu(hObject)

    handles = guidata(hObject) ;
    
    mouse = handles.mouse_menu.Value ;
    h = handles.data_holder{mouse} ;
    
    handles.session_menu.String = {'Baseline'} ;
    for i = 2 : numel(h)
        handles.session_menu.String{i} = ['Bias ' num2str(i-1)] ;
    end
%     handles.session_menu.String{end+1} = 'Bias 2+3' ;
    
end

% ========================== Callbacks ===================================

function dFF_mode_menu_Callback(hObject, eventdata, handles)
    plot_dFF(hObject)
end

function load_data_button_Callback(hObject, eventdata, handles)
    
    [data_file,data_path] = uigetfile('*.mat' , 'Choose data file' , handles.main_folder) ;
    if data_file
        
        load([data_path data_file])
        update_anlss_menu(hObject)
        update_neuron_filter_menu(hObject)
        handles.dFF_mode_menu.String = {'Mean','Trials'} ;
        handles.odor_menu.String = {'All odors','Odor 1','Odor 2','Odor 3','Odors 1 & 3'} ;
                
        handles.data_holder = d ;

        handles.session_menu.Enable = 'on' ;
        handles.neuron_menu.Enable = 'on' ;
        handles.neuron_filter_menu.Enable = 'on' ;
        handles.mouse_menu.Enable = 'on' ;
        handles.all_mice_anlss_menu.Enable = 'on' ;
        handles.one_mouse_anlss_menu.Enable = 'on' ;
        handles.dFF_mode_menu.Enable = 'on' ;
        handles.odor_menu.Enable = 'on' ;
        handles.save_hit_button.Enable = 'On' ;
        handles.save_CR_button.Enable = 'On' ;
        handles.switch_anlss_button.Enable = 'On' ;

        handles.mouse_menu.String = {} ;
        for i = 1 : numel(handles.data_holder)
            handles.mouse_menu.String{end+1} = handles.data_holder{i}(1).Parameters.mouse ;

            p = handles.data_holder{i}(end).Parameters ;
            odor_probs(1) = p.Go_odor1_prob/(p.Go_odor1_prob + p.NoGo_odor1_prob) ;
            odor_probs(2) = p.Go_odor2_prob/(p.Go_odor2_prob + p.NoGo_odor2_prob) ;
            odor_probs(3) = p.Go_odor3_prob/(p.Go_odor3_prob + p.NoGo_odor3_prob) ;
            [~,odor_order] = sort(odor_probs,'descend') ;
            odor_order(4) = 4 ;
            for j = 1 : numel(handles.data_holder{i})
                handles.data_holder{i}(j).odor_order = odor_order ;
                handles.data_holder{i}(j).organize_data()
            end

        end
        handles.mouse_menu.Value = 1 ;
        
        guidata(hObject,handles)

        mouse_menu_Callback(hObject,0,0)
        
    end
    
end

function mouse_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    
    handles.neuron_menu.Value = 1 ;
    update_session_menu(hObject)
    handles.session_menu.Value = 1 ;
    session_menu_Callback(hObject,0,0)
        
    guidata(hObject,handles)
    
end

function neuron_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;

    plot_dFF(hObject)
    
%     if handles.one_mouse_anlss_menu.UserData
%         one_mouse_anlss_menu_Callback(hObject,0,0) ;
%     end
    
end

function odor_menu_Callback(hObject, eventdata, handles)
    plot_dFF(hObject)
end

function save_anlss_ax_Callback(hObject, eventdata, handles)

    hObject.UserData = 1 ;
    if handles.all_mice_anlss_menu.UserData
        all_mice_anlss_menu_Callback(hObject, 0, 0)
    else
        one_mouse_anlss_menu_Callback(hObject, 0, 0)
    end
    hObject.UserData = 0 ;    
    
end

function save_CR_button_Callback(hObject, eventdata, handles)
    plot_single(hObject,2)
end

function save_hit_button_Callback(hObject, eventdata, handles)
    plot_single(hObject,1)
end

function session_menu_Callback(hObject, eventdata, handles)

    handles = guidata(hObject) ;
    sess = handles.session_menu.Value ;
    mouse = handles.mouse_menu.Value ;
    if sess <= numel(handles.session_menu.String)
        h = handles.data_holder{mouse}(sess) ;

        handles.first_trial.String = num2str(h.ft) ;
        handles.first_trial.Enable = 'on' ;
        handles.last_trial.String = num2str(h.lt) ;
        handles.last_trial.Enable = 'on' ;
        handles.max_trials.String = ['/' num2str(length(h.trial_types))] ;

        neuron_filter_menu_Callback(hObject)

        update_probs_panel(hObject)
    else
    end
    
%     if handles.one_mouse_anlss_menu.UserData
%         one_mouse_anlss_menu_Callback(hObject,0,0) ; 
%     else
%         all_mice_anlss_menu_Callback(hObject,0,0) ;       
%     end 
    
end

function switch_anlss_button_Callback(hObject, eventdata, handles)

    if handles.one_mouse_anlss_menu.UserData
        handles.all_mice_anlss_menu.UserData = 1 ;
        handles.one_mouse_anlss_menu.UserData = 0 ;
        all_mice_anlss_menu_Callback(hObject,0,0) ;
    else
        handles.all_mice_anlss_menu.UserData = 0 ;
        handles.one_mouse_anlss_menu.UserData = 1 ;
        one_mouse_anlss_menu_Callback(hObject,0,0) ;
    end        
        
end

% ========================== Creation functions ==========================

function all_mice_anlss_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function dFF_mode_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function mouse_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function neuron_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function neuron_filter_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function odor_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

function one_mouse_anlss_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function session_menu_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
end

