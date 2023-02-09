function X = updateResults(ModelData, System, AnalysisParam, SYS, r_full, tt)

nElem = [ModelData.NumBeamElem];
nPart = numel(ModelData);

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %Preallocate
        X.total_Aero    = zeros(6,1);
        X.total_Grav    = zeros(6,1);
        X.total_Aero_G  = zeros(6,1);
        X.total_Grav_G  = zeros(6,1);
        X.tail_Aero_G   = zeros(6,1);
        X.wing_Aero_G   = zeros(6,1);
        X.member_Aero_G = cell(nPart,1);
        X.member_Aero   = cell(nPart,1);
        X.member_Grav_G = cell(nPart,1);
        X.member_Grav   = cell(nPart,1);
        CGa             = System.CGa;
        
        %States map straight across
        X.x_f = System.x_f;
        
        %Get element orientations
        [~ , CaB] = strains2coords(ModelData, System.x_f);
        
        for jj = 1 : nPart
            
            Btot_mp =  ModelData(jj).ShapeFunctions.Integral;
            
            %THIS CODE ASSUMES THAT ALL LOADS IN 'SYS.f' ARE ATTRIBUTED TO
            %AERODYNAMIC FORCES! THIS IS WRONG!!!
            
            %'aeroForce' is in the local frame
            %X.aeroForce{jj}   = SYS.f{jj};
            %X.Local_AoA{jj}   = SYS.AoA_qs_t{jj}*180/pi;
            
            ub = cumsum(repmat(6, [1, nElem(jj)]));
            lb = [1, ub(1 : end - 1) + 1];
            
            %Convert forces from the element frame to the global frame 
            for ii = 1 : nElem(jj)
                elemInd = lb(ii) : ub(ii);
                %Element-to-global transformation matrix
                CGB_e = CGa * CaB(:,:,ModelData(jj).MatrixIndexing.NodeMidInd(ii));
                X.aeroForce_G{jj}(elemInd, 1) = [CGB_e,zeros(3);zeros(3),CGB_e]*SYS.f(elemInd);
                X.gravForce_G{jj}(elemInd, 1) = [CGB_e,zeros(3);zeros(3),CGB_e]*SYS.f_grav(elemInd);
            end
            
            %Combine forces and return in local and global frame
            X.allForce_G{jj}    = X.aeroForce_G{jj} + X.gravForce_G{jj};
            X.allForce_G_pt{jj} = Btot_mp*X.allForce_G{jj};
            X.member_Aero{jj}   =                              SYS.Btot_r*SYS.f;
            X.member_Aero_G{jj} = [CGa zeros(3);zeros(3) CGa]* SYS.Btot_r*SYS.f;
            X.member_Grav{jj}   =                              SYS.Btot_r*SYS.f_grav + SYS.F_r*SYS.f_grav_pt ;
            X.member_Grav_G{jj} = [CGa zeros(3);zeros(3) CGa]*(SYS.Btot_r*SYS.f_grav + SYS.F_r*SYS.f_grav_pt);
            
            X.total_Aero   = X.total_Aero   + SYS.Btot_r*SYS.f;
            X.total_Grav   = X.total_Grav   + SYS.Btot_r*SYS.f_grav + SYS.F_r*SYS.f_grav_pt;
            X.total_Aero_G = X.total_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r*SYS.f;
            X.total_Grav_G = X.total_Grav_G + [CGa zeros(3);zeros(3) CGa]*(SYS.Btot_r*SYS.f_grav + SYS.F_r*SYS.f_grav_pt);
            %TODO - Remove hardcoded terms for the wing/tail/etc.
            if jj == 1 || jj == 2
                X.wing_Aero_G = X.wing_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r*SYS.f;
            elseif jj == 3 || jj == 4 || jj == 5
                X.tail_Aero_G = X.tail_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r*SYS.f;
            end
            X.res = norm(r_full);
        end
%         if Sim.rb_flag == 1
%             X.total_Grav   = X.total_Grav   +                              Matrices.M_rb*[CGa' zeros(3); zeros(3) CGa']*Sim.grav_vec;
%             X.total_Grav_G = X.total_Grav_G + [CGa zeros(3); zeros(3) CGa]*Matrices.M_rb*[CGa' zeros(3); zeros(3) CGa']*Sim.grav_vec;
%         end
        
    otherwise        
        
        if Sim.Soln == 2 || Sim.Soln == 2.1 %Dynamic
            X.total_Aero(:,tt+1)   = zeros(6,1);
            X.total_Aero_G(:,tt+1) = zeros(6,1);
            X.tail_Aero_G(:,tt+1)  = zeros(6,1);
            X.wing_Aero_G(:,tt+1)  = zeros(6,1);
            
            CGa                    = Quat2Rot(X_p1.x_qa_p1);
            
            for jj = 1:nPart
                X.x_f{jj}(:,tt+1)     = X_p1.x_f_p1{jj};
                X.x_v{jj}(:,tt+1)     = X_p1.x_v_p1{jj};
                X.x_q{jj}(:,tt+1)     = X_p1.x_q_p1{jj};
                X.x_p{jj}(:,tt+1)     = X_p1.x_p_p1{jj};
                X.x_x{jj}(:,tt+1)     = X_p1.x_x_p1{jj};
                
                X.x_f_dot{jj}(:,tt+1) = X_p1.x_f_dot_p1{jj};
                X.x_v_dot{jj}(:,tt+1) = X_p1.x_v_dot_p1{jj};
                X.x_q_dot{jj}(:,tt+1) = X_p1.x_q_dot_p1{jj};
                X.x_p_dot{jj}(:,tt+1) = X_p1.x_p_dot_p1{jj};
                X.x_x_dot{jj}(:,tt+1) = X_p1.x_x_dot_p1{jj};
                
                %                 KE(tt+1) = KE(tt+1) + x_v{jj}(:,tt+1)' * M{jj}  * x_v{jj}(:,tt+1)/2;
                %                 PE(tt+1) = PE(tt+1) + x_f{jj}(:,tt+1)' * T1{jj} * x_f{jj}(:,tt+1)/2;
                %
                X.Force_Store{jj}(:,tt+1) = SYS.f{jj};
                if Sim.aero_flag
                    X.Local_AoA{jj}(:,tt+1)   = SYS.AoA_qs_t{jj}*180/pi;
                    
                    X.total_Aero(:,tt+1)   = X.total_Aero(:,tt+1)     + SYS.Btot_r{jj}*SYS.f{jj};
                    X.total_Aero_G(:,tt+1) = X.total_Aero_G(:,tt+1)   + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
                    if jj == 1 || jj == 2
                        X.wing_Aero_G(:,tt+1) = X.wing_Aero_G(:,tt+1) + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
                    elseif jj == 3 || jj == 4
                        X.tail_Aero_G(:,tt+1) = X.tail_Aero_G(:,tt+1) + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
                    end
                end
            end
            
            if iscell(r_full)
                norm_r_full = [];
                for ii = 1:length(r_full)
                    norm_r_full(ii) = norm(r_full{ii});
                end
                norm_r_full = norm(norm_r_full);
            else
                norm_r_full = norm(r_full);
            end
            
            X.residual(tt+1) = norm_r_full;
            
            if Sim.rb_flag
                X.x_va(:,tt+1)     = X_p1.x_va_p1;
                X.x_qa(:,tt+1)     = X_p1.x_qa_p1;
                X.x_pa(:,tt+1)     = X_p1.x_pa_p1;
                
                X.x_vg(:,tt+1)     = [CGa zeros(3);zeros(3) CGa]*X_p1.x_va_p1;
                
                X.RBEulerAngles(:,tt+1) = Quat2Euler(X_p1.x_qa_p1);%SpinCalc('QtoEA213',X_p1.x_qa_p1',[],1);
                
                X.x_va_dot(:,tt+1) = X_p1.x_va_dot_p1;
                X.x_qa_dot(:,tt+1) = X_p1.x_qa_dot_p1;
                X.x_pa_dot(:,tt+1) = X_p1.x_pa_dot_p1;
            end
            
            return
            
        end
        
        if Sim.Soln == 2.2
            X.x(:,tt+1)         = X_p1.x_p1;
            X.x_dot(:,tt+1)     = X_p1.x_dot_p1;
        end
        
end

end