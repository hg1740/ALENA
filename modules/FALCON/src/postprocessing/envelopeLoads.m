function envelopeLoads(ListToEnvelope,Listof1g,ListofLoadsShift,LoadFactor,Member,SaveAsFile,LinFlag,ListToEnvelope90)

fprintf(['\nPost-processing ',num2str(length(ListToEnvelope)),' simulations...\n'])
if LinFlag
    fprintf('\nEnveloping Linear Results...\n')
else
    fprintf('\nEnveloping Nonlinear Results...\n')
end

RTC_counter = 0;

for ff = 1:length(ListToEnvelope)
    fprintf(['\nProcessing gust result ',num2str(ff),' of ',num2str(length(ListToEnvelope))])
    
    load(ListToEnvelope{ff})
    if isempty(Listof1g{ff})
        Loads1g = zeros(size(x_out1.x_f{Member}));
    else
        x1g = load(Listof1g{ff});
        Loads1g = repmat(x1g.x_out.x_f{Member},1,size(x_out1.x_f{Member},2));
    end
    if isempty(ListofLoadsShift{ff})
        LoadsShift = zeros(size(x_out1.x_f{Member}));
    else
        xLS        = load(ListofLoadsShift{ff});
        LoadsShift = xLS.x_out1.x_f{Member};
    end

    if ff == 1
        envelope.env1D.Loads.Max    = ones( size(x_out1.x_f{Member},1),1)*-1e99;
        envelope.env1D.Loads.Min    = ones( size(x_out1.x_f{Member},1),1)* 1e99;
        envelope.env1D.Loads.MaxInc = ones( size(x_out1.x_f{Member},1),1)*-1e99;
        envelope.env1D.Loads.MinInc = ones( size(x_out1.x_f{Member},1),1)* 1e99;
        envelope.env1D.Case.Max     = cell( size(x_out1.x_f{Member},1),1);
        envelope.env1D.Case.Min     = cell( size(x_out1.x_f{Member},1),1);
        envelope.env1D.Time.Max     = ones( size(x_out1.x_f{Member},1),1)* NaN;
        envelope.env1D.Time.Min     = ones( size(x_out1.x_f{Member},1),1)* NaN;
%         envelope.env2D.Loads.theta  = [];
%         envelope.env2D.Loads.MxvMy  = [];
        envelope.env2D.Loads.Mx     = [];
        envelope.env2D.Loads.My     = [];
    end
    
    if LinFlag
        loads     = [Loads1g + LoadFactor*(x_out1.x_f{Member} - LoadsShift) , Loads1g - LoadFactor*(x_out1.x_f{Member} - LoadsShift)];
    else
        loads     =  Loads1g + LoadFactor*(x_out1.x_f{Member} - LoadsShift);
    end
    
%     envelope.env2D.Loads.MxvMy = [envelope.env2D.Loads.MxvMy (loads(4,:).^2+loads(5,:).^2).^0.5];
%     envelope.env2D.Loads.theta = [envelope.env2D.Loads.theta atan((loads(5,:)-loads(5,1))./(loads(4,:)-loads(4,1)))      ];
    envelope.env2D.Loads.Mx    = [envelope.env2D.Loads.Mx loads(4,:)];
    envelope.env2D.Loads.My    = [envelope.env2D.Loads.My loads(5,:)];
        
    [max_val,max_ind] = max([envelope.env1D.Loads.Max , loads],[],2);
    [min_val,min_ind] = min([envelope.env1D.Loads.Min , loads],[],2);
    
    envelope.env1D.PerCase.Max(:,ff) = max(loads,[],2);
    envelope.env1D.PerCase.Min(:,ff) = min(loads,[],2);
    if LinFlag
        envelope.env1D.PerCase.MaxInc(:,ff) = envelope.env1D.PerCase.Max(:,ff) - Loads1g(:,1);
        envelope.env1D.PerCase.MinInc(:,ff) = envelope.env1D.PerCase.Min(:,ff) - Loads1g(:,1);
    else
        envelope.env1D.PerCase.MaxInc(:,ff) = envelope.env1D.PerCase.Max(:,ff) - x_out1.x_f{Member}(:,1);
        envelope.env1D.PerCase.MinInc(:,ff) = envelope.env1D.PerCase.Min(:,ff) - x_out1.x_f{Member}(:,1);
    end
    envelope.env1D.PerCase.CaseName{1,ff} = ListToEnvelope{ff};
    
    for pp = 1:length(envelope.env1D.Loads.Max)
        if max_ind(pp) > 1
            indToUse = max_ind(pp) - 1;
            if indToUse > size(x_out1.x_f{Member},2)
                indToUse = indToUse - size(x_out1.x_f{Member},2);
                envelope.env1D.Case.Max{pp}     = [ListToEnvelope{ff} '(Negative)'];
            else
                envelope.env1D.Case.Max{pp}     =  ListToEnvelope{ff};
            end
            
            envelope.env1D.Loads.Max(pp)    = max_val(pp);
            if LinFlag
                envelope.env1D.Loads.MaxInc(pp) = max_val(pp) - Loads1g(pp,1);
            else
                envelope.env1D.Loads.MaxInc(pp) = max_val(pp) - x_out1.x_f{Member}(pp,1);
            end
            envelope.env1D.Time.Max(pp)     = indToUse;
        end
        if min_ind(pp) > 1
            indToUse = min_ind(pp) - 1;
            if indToUse > size(x_out1.x_f{Member},2)
                indToUse = indToUse - size(x_out1.x_f{Member},2);
                envelope.env1D.Case.Min{pp}     = [ListToEnvelope{ff} '(Negative)'];
            else
                envelope.env1D.Case.Min{pp}     =  ListToEnvelope{ff};
            end
            
            envelope.env1D.Loads.Min(pp)    = min_val(pp);
            if LinFlag
                envelope.env1D.Loads.MinInc(pp) = min_val(pp) - Loads1g(pp,1);
            else
                envelope.env1D.Loads.MinInc(pp) = min_val(pp) - x_out1.x_f{Member}(pp,1);
            end
            envelope.env1D.Time.Min(pp)     = indToUse;
        end
    end
    if LinFlag && ~isempty(ListToEnvelope90{ff})
        fprintf(['\nProcessing gust result ',num2str(ff),' of ',num2str(length(ListToEnvelope)),' (RTC gusts)'])
        
        RTC_counter = RTC_counter + 1;
        
        LoadsVert = load(ListToEnvelope{ff} );
        LoadsLat  = load(ListToEnvelope90{ff});
        if isempty(Listof1g{ff})
            Loads1g = zeros(size(LoadsVert.x_out1.x_f{Member}));
        else
            x1g = load(Listof1g{ff});
            Loads1g = repmat(x1g.x_out.x_f{Member},1,size(LoadsVert.x_out1.x_f{Member},2));
        end
        if isempty(ListofLoadsShift{ff})
            LoadsShift = zeros(size(LoadsVert.x_out1.x_f{Member}));
        else
            xLS        = load(ListofLoadsShift{ff});
            LoadsShift = xLS.x_out1.x_f{Member};
        end
        
        if RTC_counter == 1
            envelope.RTC.Loads.Max    = ones(size(LoadsVert.x_out1.x_f{Member},1),1)*-1e99;
            envelope.RTC.Loads.Min    = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* 1e99;
            envelope.RTC.Loads.MaxInc = ones(size(LoadsVert.x_out1.x_f{Member},1),1)*-1e99;
            envelope.RTC.Loads.MinInc = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* 1e99;
            envelope.RTC.Theta.Max    = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* NaN;
            envelope.RTC.Theta.Min    = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* NaN;
            envelope.RTC.Time.Max     = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* NaN;
            envelope.RTC.Time.Min     = ones(size(LoadsVert.x_out1.x_f{Member},1),1)* NaN;
            envelope.RTC.Case.Max     = cell(size(LoadsVert.x_out1.x_f{Member},1),1);
            envelope.RTC.Case.Min     = cell(size(LoadsVert.x_out1.x_f{Member},1),1);
        end
        
        loads0    = LoadFactor*(LoadsVert.x_out1.x_f{Member} - LoadsShift);
        loads90   = LoadFactor*(LoadsLat.x_out1.x_f{Member}  - LoadsShift);
        
        theta = atan(loads90./loads0)*180/pi;
        theta(isnan(theta)) = 0;
        loadsRTCMax  = Loads1g + (loads0.^2+loads90.^2).^0.5;
        loadsRTCMin  = Loads1g - (loads0.^2+loads90.^2).^0.5;
        
        [max_val,max_ind] = max([envelope.RTC.Loads.Max , loadsRTCMax, loadsRTCMin],[],2);
        [min_val,min_ind] = min([envelope.RTC.Loads.Min , loadsRTCMax, loadsRTCMin],[],2);
        
        envelope.RTC.PerCase.Max(:,RTC_counter)      = max([loadsRTCMax, loadsRTCMin],[],2);
        envelope.RTC.PerCase.Min(:,RTC_counter)      = min([loadsRTCMax, loadsRTCMin],[],2);
        
        envelope.RTC.PerCase.MaxInc(:,RTC_counter)   = envelope.RTC.PerCase.Max(:,RTC_counter) - Loads1g(:,1);
        envelope.RTC.PerCase.MinInc(:,RTC_counter)   = envelope.RTC.PerCase.Min(:,RTC_counter) - Loads1g(:,1);

        envelope.RTC.PerCase.CaseName{1,RTC_counter} = ListToEnvelope{ff};
        
        for pp = 1:length(envelope.RTC.Loads.Max)
            if max_ind(pp) > 1
                indToUse = max_ind(pp) - 1;
                if indToUse > size(LoadsVert.x_out1.x_f{Member},2)
                    indToUse = indToUse - size(LoadsVert.x_out1.x_f{Member},2);
                    envelope.RTC.Theta.Max(pp)    = -theta(pp,max_ind(pp)-1);
                    envelope.RTC.Case.Max{pp}     = [ListToEnvelope{ff} '(Negative)'];
                else
                    envelope.RTC.Theta.Max(pp)   =   theta(pp,max_ind(pp)-1);
                    envelope.RTC.Case.Max{pp}     =  ListToEnvelope{ff};
                end
                
                envelope.RTC.Loads.Max(pp)    = max_val(pp);
                envelope.RTC.Loads.MaxInc(pp) = max_val(pp) - Loads1g(pp,1);
                envelope.RTC.Time.Max(pp)     = max_ind(pp)-1;
            end
            if min_ind(pp) > 1
                indToUse = min_ind(pp) - 1;
                if indToUse > size(LoadsVert.x_out1.x_f{Member},2)
                    indToUse = indToUse - size(LoadsVert.x_out1.x_f{Member},2);
                    envelope.RTC.Theta.Min(pp)    = -theta(pp,indToUse);
                    envelope.RTC.Case.Min{pp}     = [ListToEnvelope{ff} '(Negative)'];
                else
                    envelope.RTC.Theta.Min(pp)   =   theta(pp,indToUse);
                    envelope.RTC.Case.Max{pp}     =  ListToEnvelope{ff};
                end
                
                envelope.RTC.Loads.Min(pp)    = min_val(pp);
                envelope.RTC.Loads.MinInc(pp) = min_val(pp) - Loads1g(pp,1);
                envelope.RTC.Time.Min(pp)     = min_ind(pp)-1;
            end
        end
    end
end

fprintf('\nCalculating 2D envelopes...\n')

k = convhull(envelope.env2D.Loads.Mx,envelope.env2D.Loads.My);
envelope.env2D.Loads.Mx    = envelope.env2D.Loads.Mx(k);
envelope.env2D.Loads.My    = envelope.env2D.Loads.My(k);
envelope.env2D.Loads.aveMx = mean(envelope.env2D.Loads.Mx);
envelope.env2D.Loads.aveMy = mean(envelope.env2D.Loads.My);
envelope.env2D.Loads.MxvMy = ((envelope.env2D.Loads.Mx-envelope.env2D.Loads.aveMx).^2 + (envelope.env2D.Loads.My-envelope.env2D.Loads.aveMy).^2).^0.5;
envelope.env2D.Loads.theta = atan2((envelope.env2D.Loads.My-envelope.env2D.Loads.aveMy),(envelope.env2D.Loads.Mx-envelope.env2D.Loads.aveMx));
[envelope.env2D.Loads.theta,ind] = sort(envelope.env2D.Loads.theta);
envelope.env2D.Loads.MxvMy = envelope.env2D.Loads.MxvMy(ind);

Mx_reorder   = envelope.env2D.Loads.Mx(ind)-envelope.env2D.Loads.aveMx; 
My_reorder   = envelope.env2D.Loads.My(ind)-envelope.env2D.Loads.aveMy; 
Mx_crossover = [Mx_reorder(1) Mx_reorder(end)];
My_crossover = [My_reorder(1) My_reorder(end)];

MxvMy_interp = interp1(My_crossover,Mx_crossover,0);

envelope.env2D.Loads.theta = [-pi               envelope.env2D.Loads.theta pi               ];
envelope.env2D.Loads.MxvMy = [abs(MxvMy_interp) envelope.env2D.Loads.MxvMy abs(MxvMy_interp)];

theta_refined = -pi:2*pi/10000:pi;

[~,ind] = unique(envelope.env2D.Loads.theta);
Mx_refined = interp1(envelope.env2D.Loads.theta(ind),envelope.env2D.Loads.MxvMy(ind).*cos(envelope.env2D.Loads.theta(ind)),theta_refined,'linear');
My_refined = interp1(envelope.env2D.Loads.theta(ind),envelope.env2D.Loads.MxvMy(ind).*sin(envelope.env2D.Loads.theta(ind)),theta_refined,'linear');
envelope.env2D.Loads.MxvMy = (Mx_refined.^2+My_refined.^2).^0.5;
theta_refined_check = atan2(My_refined,Mx_refined);
theta_refined_check(end) = abs(theta_refined_check(end));
% plot(Mx_refined,My_refined,'k')
% hold on
for ii = 1:length(theta_refined)
    if ~isempty(find(theta_refined(ii)==theta_refined_check,1))
        Mx_refined_2(ii) = Mx_refined(theta_refined(ii)==theta_refined_check);
        My_refined_2(ii) = My_refined(theta_refined(ii)==theta_refined_check);
    else
        [theta_p] = min(theta_refined_check(theta_refined_check>theta_refined(ii)));
        [theta_m] = max(theta_refined_check(theta_refined_check<theta_refined(ii)));
        ind_p     = find(theta_refined_check == theta_p);
        ind_m     = find(theta_refined_check == theta_m);
        Mx_p = Mx_refined(ind_p);
        My_p = My_refined(ind_p);
        Mx_m = Mx_refined(ind_m);
        My_m = My_refined(ind_m);
        coeffs = [Mx_p 1; Mx_m 1]\[My_p; My_m];
        a = coeffs(1); b = coeffs(2);
        Mx_refined_2(ii) =   b/(tan(theta_refined(ii))-a)    ;
        My_refined_2(ii) = a*b/(tan(theta_refined(ii))-a) + b;
%         scatter(Mx_p,My_p,'b.')
%         scatter(Mx_m,My_m,'b.')
    end
%     scatter(Mx_refined_2(ii),My_refined_2(ii),'r.')
end
% plot(Mx_refined_2,My_refined_2,'r--')

envelope.env2D.Loads.theta = theta_refined;
envelope.env2D.Loads.MxvMy = (Mx_refined_2.^2+My_refined_2.^2).^0.5;

fprintf('\nSaving envelope...\n')
save(SaveAsFile,'envelope')
fprintf(['\nEnvelope Saved in ''',strrep(SaveAsFile,'\','\\'),'''!\n'])