function IND = getSystemInd(Matrices, Sim, AeroFlag, AnalysisParam)
%getSystemInd Returns the linear index numbers for each component in the
%system.

%How many elements in each component?
nElem = Matrices.n_elem;

%How many components?
nComp = numel(nElem);

%Indexing bounds
ub = cumsum(nElem * 6);
lb = [1, ub(1 : end - 1) + 1];

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %Linear index numbers for each component
        IND.x_f_ind = arrayfun(@(i) lb(i) : ub(i), 1 : nComp, 'Unif', false);

    otherwise
        
        if Sim.Soln == 2
            if Sim.aero_flag && Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_p_ind{jj}(end);
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_x_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_p_ind{jj}(end);
                    end
                    IND.flex_ind{jj}    = IND.x_v_ind{jj}(1):IND.x_x_ind{jj}(end);
                    IND.BC1_tot{jj}     = blkdiag(Matrices.R1_bc{jj},Matrices.R2_bc{jj},eye(  Matrices.n_elem(jj)*4),eye(  Matrices.n_elem(jj)*3),eye(  Matrices.n_elem(jj)*2));
                    IND.BC2_tot{jj}     = blkdiag(Matrices.R3_bc{jj},Matrices.R4_bc{jj},zeros(Matrices.n_elem(jj)*4),zeros(Matrices.n_elem(jj)*3),zeros(Matrices.n_elem(jj)*2));
                end
                IND.x_va_ind    = 1:6;
                IND.x_qa_ind    = 1:4;
                IND.x_pa_ind    = 1:3;
                IND.x_va_ind = IND.x_va_ind + IND.x_x_ind{length(Matrices.n_elem)}(end);
                IND.x_qa_ind = IND.x_qa_ind + IND.x_va_ind(end);
                IND.x_pa_ind = IND.x_pa_ind + IND.x_qa_ind(end);
            elseif  Sim.aero_flag && ~Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_p_ind{jj}(end);
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_x_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_p_ind{jj}(end);
                    end
                    IND.flex_ind{jj}    = IND.x_v_ind{jj}(1):IND.x_x_ind{jj}(end);
                    IND.BC1_tot{jj}     = blkdiag(Matrices.R1_bc{jj},Matrices.R2_bc{jj},eye(  Matrices.n_elem(jj)*4),eye(  Matrices.n_elem(jj)*3),eye(  Matrices.n_elem(jj)*2));
                    IND.BC2_tot{jj}     = blkdiag(Matrices.R3_bc{jj},Matrices.R4_bc{jj},zeros(Matrices.n_elem(jj)*4),zeros(Matrices.n_elem(jj)*3),zeros(Matrices.n_elem(jj)*2));
                end
                IND.x_va_ind    = [];
                IND.x_qa_ind    = [];
                IND.x_pa_ind    = [];
            elseif ~Sim.aero_flag &&  Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_p_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    end
                    IND.flex_ind{jj}    = IND.x_v_ind{jj}(1):IND.x_p_ind{jj}(end);
                    IND.BC1_tot{jj}     = blkdiag(Matrices.R1_bc{jj},Matrices.R2_bc{jj},eye(  Matrices.n_elem(jj)*4),eye(  Matrices.n_elem(jj)*3));
                    IND.BC2_tot{jj}     = blkdiag(Matrices.R3_bc{jj},Matrices.R4_bc{jj},zeros(Matrices.n_elem(jj)*4),zeros(Matrices.n_elem(jj)*3));
                end
                IND.x_va_ind    = 1:6;
                IND.x_qa_ind    = 1:4;
                IND.x_pa_ind    = 1:3;
                IND.x_va_ind = IND.x_va_ind + IND.x_p_ind{length(Matrices.n_elem)}(end);
                IND.x_qa_ind = IND.x_qa_ind + IND.x_va_ind(end);
                IND.x_pa_ind = IND.x_pa_ind + IND.x_qa_ind(end);
            elseif ~Sim.aero_flag && ~Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_p_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    end
                    IND.flex_ind{jj}    = IND.x_v_ind{jj}(1):IND.x_p_ind{jj}(end);
                    IND.BC1_tot{jj}     = blkdiag(Matrices.R1_bc{jj},Matrices.R2_bc{jj},eye(  Matrices.n_elem(jj)*4),eye(  Matrices.n_elem(jj)*3));
                    IND.BC2_tot{jj}     = blkdiag(Matrices.R3_bc{jj},Matrices.R4_bc{jj},zeros(Matrices.n_elem(jj)*4),zeros(Matrices.n_elem(jj)*3));
                end
                IND.x_va_ind    = [];
                IND.x_qa_ind    = [];
                IND.x_pa_ind    = [];
            end
        elseif Sim.Soln == 2.1
            if Sim.aero_flag && Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6                          ;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end)  ;
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4)                         ;
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end)  ;
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_f_ind{jj}(end)  ;
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_x_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end)  ;
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_p_ind{jj-1}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end)  ;
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_f_ind{jj}(end)  ;
                    end
                end
                IND.x_va_ind    = 1:6;
                IND.x_qa_ind    = 1:4;
                IND.x_pa_ind    = 1:3;
                IND.x_va_ind = IND.x_va_ind + IND.x_x_ind{Matrices.n_elem(end)}(end);
                IND.x_qa_ind = IND.x_qa_ind + IND.x_p_ind{Matrices.n_elem(end)}(end);
                IND.x_pa_ind = IND.x_pa_ind + IND.x_qa_ind(end);
            elseif  Sim.aero_flag && ~Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6                          ;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end)  ;
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4)                         ;
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end)  ;
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_f_ind{jj}(end)  ;
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_x_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end)  ;
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_p_ind{jj-1}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end)  ;
                        IND.x_x_ind{jj} = (1:Matrices.n_elem(jj)*2) + IND.x_f_ind{jj}(end)  ;
                    end
                end
                IND.x_va_ind    = [];
                IND.x_qa_ind    = [];
                IND.x_pa_ind    = [];
            elseif ~Sim.aero_flag &&  Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_p_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    end
                end
                IND.x_va_ind    = 1:6;
                IND.x_qa_ind    = 1:4;
                IND.x_pa_ind    = 1:3;
                IND.x_va_ind = IND.x_va_ind + IND.x_p_ind{length(Matrices.n_elem)}(end);
                IND.x_qa_ind = IND.x_qa_ind + IND.x_va_ind(end);
                IND.x_pa_ind = IND.x_pa_ind + IND.x_qa_ind(end);
            elseif ~Sim.aero_flag && ~Sim.rb_flag
                for jj = 1:length(Matrices.n_elem)
                    if jj == 1
                        IND.x_v_ind{jj} =  1:Matrices.n_elem(jj)*6;
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    else
                        IND.x_v_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_p_ind{jj-1}(end);
                        IND.x_f_ind{jj} = (1:Matrices.n_elem(jj)*6) + IND.x_v_ind{jj}(end);
                        IND.x_q_ind{jj} = (1:Matrices.n_elem(jj)*4) + IND.x_f_ind{jj}(end);
                        IND.x_p_ind{jj} = (1:Matrices.n_elem(jj)*3) + IND.x_q_ind{jj}(end);
                        IND.x_x_ind{jj} = [];
                    end
                end
                IND.x_va_ind    = [];
                IND.x_qa_ind    = [];
                IND.x_pa_ind    = [];
            end
        elseif Sim.Soln == 2.2
            IND.x_v_ind{1} =  1:Matrices.n_elem(1)*6;
            IND.x_f_ind{1} = (1:Matrices.n_elem(1)*6) + IND.x_v_ind{1}(end);
        end
        
end

end