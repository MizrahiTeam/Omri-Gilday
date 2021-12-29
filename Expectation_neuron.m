classdef Expectation_neuron < handle & matlab.mixin.Heterogeneous
    
    % ============ Fields =================================
    
    properties
        
        Parameters
        
        dFF
        F0
        F0_s
        resps
        mean_resps
        valid_trials
        valid_trials_clssf
        invalid = 0 ;
        responsiveness
        responsive = 0 ;
        odor_responsiveness;
        pref_sound % 1 for Go 0 for NoGo -1 for none-selective
        SPI
        exp_DI = zeros(1,3) ;
        exp_DI_diff = 0 ;
        
        odor_resps  
        inhibited
       
    end
    
    properties (Constant = true)        
        
        frames_before_stim = 16 ;   % Number of frames to be taken before stimulus
        frames_after_stim = 37  ;   % Number of frames to be taken after stimulus
        baseline_frames = 7     ;   % Number of frames to be taken as baseline (starting from first frame) 
        frames_for_amp = 7     ;    % Number of frames to be taken into account for ampltude calculation
        MAX_F0_s = 0.2 ;
        MIN_AMP_SIG = 0.05 ;
        MIN_TRIALS = 5 ;
        SPI_thresh = 0.1 ;
                
    end
    
    % ============ Methods ================================
    
    methods
        
        function obj = Expectation_neuron(BOT,trial_types,trial_odors,tSOF,tSOS,ft_all,lt_all,Parameters)
            
            obj.Parameters = Parameters ;
            obj.make_dFF(BOT,trial_types,trial_odors,tSOF,tSOS)
            obj.find_valid_trials(ft_all,lt_all)
            obj.check_responsiveness() ;
            obj.sound_preference() ;
            obj.expectation_DI() ;
            
        end
        
        function make_dFF(obj,BOT,trial_types,trial_odors,tSOF,tSOS)
            
            n_frames = obj.frames_before_stim + obj.frames_after_stim ;
            base_st = obj.frames_before_stim - obj.baseline_frames + 1 ; % first frame for odor baseline
            base_en = obj.frames_before_stim ; % last frame for odor baseline
            sound_st = obj.frames_before_stim + 1 ; % first frame for stim response
            sound_en = obj.frames_before_stim + obj.frames_for_amp ; % last frame for stim response
            odor_st = obj.frames_before_stim - floor(7.2*obj.Parameters.odor_dur*1e-3) ;
            odor_basline_st = obj.frames_before_stim - 2*floor(7.2*obj.Parameters.odor_dur*1e-3) ;
            
            for odor = 1 : 3
                for type = 1 : 4
                    n_trials = sum(trial_odors==odor & trial_types==type) ;
                    obj.dFF{odor,type} = zeros(n_trials,n_frames) ;
                    obj.F0{odor,type} = zeros(1,n_trials) ;
                    obj.F0_s{odor,type} = zeros(1,n_trials) ;
                    obj.resps{odor,type} = zeros(1,n_trials) ;
                    obj.odor_resps{odor,type} = zeros(1,n_trials) ;
                end
            end
            
            counts = ones(3,4) ;

            for trial = 1 : length(trial_types)
                if trial_types(trial) < 7
                    type = trial_types(trial) ;
                    odor = trial_odors(trial) ;
                    x = counts(odor,type) ;
                    be = find(tSOF<tSOS(trial),obj.frames_before_stim,'last') ;
                    af = find(tSOF>tSOS(trial),obj.frames_after_stim,'first') ;
                    trace = BOT([be ; af]) ;
                    obj.F0{odor,type}(x) = mean(trace(base_st:base_en)) ;
                    obj.F0_s{odor,type}(x) = std(trace(base_st:base_en))./ obj.F0{odor,type}(x) ;
                    y = trace./(obj.F0{odor,type}(x)*ones(n_frames,1)) - 1 ;
                    obj.dFF{odor,type}(x,:) = y ;
                    obj.resps{odor,type}(x) = mean(y(sound_st:sound_en)) ;
                    counts(odor,type) = x + 1 ;
                    odor_F0 = mean(trace(odor_basline_st:(odor_st-1))) ;
                    y_odor = trace./(odor_F0*ones(n_frames,1)) - 1 ;
                    obj.odor_resps{odor,type}(x) = mean(y_odor(odor_st:(sound_st-1))) ;
                end
            end
            
            for odor = 1 : 3
                for type = 1 : 4
                    obj.mean_resps(odor,type) = mean(obj.resps{odor,type}) ;
                end
            end
            
        end
        
        function find_valid_trials(obj,ft_all,lt_all)
            
            obj.valid_trials = cell(3,4) ;
            obj.valid_trials_clssf = cell(3,4) ;
            n_valid_trials = zeros(3,4) ;
                        
            for odor = 1 : 3
                for type = 1 : 4
                    obj.valid_trials{odor,type} = zeros(size(obj.F0_s{odor,type})) ;
                    obj.valid_trials{odor,type}(ft_all(odor,type):lt_all(odor,type)) = 1 ;
                    obj.valid_trials_clssf{odor,type} = obj.valid_trials{odor,type} ;
                    obj.valid_trials{odor,type}(obj.F0_s{odor,type}>obj.MAX_F0_s) = 0 ;
                    obj.valid_trials{odor,type} = logical(obj.valid_trials{odor,type}) ;
                    obj.valid_trials_clssf{odor,type} = logical(obj.valid_trials_clssf{odor,type}) ;
                    n_valid_trials(odor,type) = sum(obj.valid_trials{odor,type}) ;
                end
            end
            
            if sum(sum(n_valid_trials(:,1:2)<obj.MIN_TRIALS))
                obj.invalid = 1 ;
            end
            
        end
        
        function check_responsiveness(obj)
            
            obj.responsiveness = zeros(3,2) ;
            obj.odor_responsiveness = zeros(3,2) ;
            
            sound_st = obj.frames_before_stim + 1 ; % first frame for stim response
            sound_en = obj.frames_before_stim + obj.frames_for_amp ; % last frame for stim response
            
            if ~obj.invalid
                for type = 1:2
                    for odor = 1 : 3
                        try
                            p = signrank(obj.get_responses(odor,type),0,'Tail','right') ;
                            mean_trace = mean(obj.dFF{odor,type},1) ;
                            amp = max(mean_trace(sound_st:sound_en)) ;
                            if (p<0.05) && amp >= obj.MIN_AMP_SIG
                                obj.responsiveness(odor,type) = 1 ;
                            end
                            p = signrank(obj.get_responses(odor,type),0,'Tail','left') ;
                            amp = min(mean_trace(sound_st:sound_en)) ;
                            if (p<0.05) && amp <= -obj.MIN_AMP_SIG
                                obj.invalid = 1 ;
                                obj.inhibited = 1 ;
                            end
                        catch
                        end
                    end
                end
            end
            
            if sum(sum(obj.responsiveness(:,1:2)))
                obj.responsive = 1 ;
            end
            
        end
        
        function sound_preference(obj)

            if ~obj.invalid && obj.responsive
                
                Hit = max([mean(obj.get_responses(1,1)) mean(obj.get_responses(2,1)) mean(obj.get_responses(3,1))]) ;
                CR = max([mean(obj.get_responses(1,2)) mean(obj.get_responses(2,2)) mean(obj.get_responses(3,2))]) ;
                if Hit > CR
                    obj.pref_sound = 1 ;
                elseif Hit < CR
                    obj.pref_sound = 0 ;
                else
                    obj.pref_sound = -1 ;
                end

            end
        
        end

        function expectation_DI(obj) 
            
            if ~obj.invalid && obj.responsive
                exp_hit = obj.get_responses(1,1) ;
                exp_CR = obj.get_responses(3,2) ;
                neut_hit = obj.get_responses(2,1) ;
                neut_CR = obj.get_responses(2,2) ;
                unexp_hit = obj.get_responses(3,1) ;
                unexp_CR = obj.get_responses(1,2) ;

                
                [~,~,~,temp(1)] = perfcurve([ones(size(exp_hit)) zeros(size(exp_CR))],[exp_hit exp_CR],1) ;
                

                [~,~,~,temp(2)] = perfcurve([ones(size(neut_hit)) zeros(size(neut_CR))],[neut_hit neut_CR],1) ;

                [~,~,~,temp(3)] = perfcurve([ones(size(unexp_hit)) zeros(size(unexp_CR))],[unexp_hit unexp_CR],1) ;
                
                for pair = 1 : 3
                    
                    obj.exp_DI(pair) = 0.5 + abs(0.5-temp(pair)) ;
                    
                end

                obj.exp_DI_diff = obj.exp_DI(3) - obj.exp_DI(1) ;
            end
            
        end
        
        function spec_resps = get_responses(obj,odor,type,is_clssf)
            
            if odor<=3
                if nargin > 3 && strcmp(is_clssf,'classifier')
                    spec_resps = obj.resps{odor,type}(obj.valid_trials_clssf{odor,type}) ;
                else
                    spec_resps = obj.resps{odor,type}(obj.valid_trials{odor,type}) ;
                end
            else
                spec_resps = [] ;
                for odor = 1 : 3
                    spec_resps = [spec_resps obj.resps{odor,type}(obj.valid_trials{odor,type})] ;
                end
            end
            
        end
                
        function valid_dFF = get_dFF(obj,odor,type)
            
            valid_dFF = obj.dFF{odor,type}(obj.valid_trials{odor,type},:) ;
            
        end
                        
    end
    
end