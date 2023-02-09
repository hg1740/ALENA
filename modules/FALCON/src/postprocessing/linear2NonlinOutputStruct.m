function x_out = linear2NonlinOutputStruct(x,Matrices,Sim)

Sim.Ind      = getSystemInd(Matrices,Sim,[]);

for jj = 1:length(Matrices.n_elem)
    x_out.x_f{jj} = x(:,Sim.Ind.x_f_ind{jj})';
    x_out.x_v{jj} = x(:,Sim.Ind.x_v_ind{jj})';
    x_out.x_p{jj} = x(:,Sim.Ind.x_p_ind{jj})';
    x_out.x_q{jj} = x(:,Sim.Ind.x_q_ind{jj})';
    if Sim.aero_flag
        x_out.x_x{jj} = x(:,Sim.Ind.x_x_ind{jj})';
    end
end
if Sim.rb_flag
    x_out.x_va = x(:,Sim.Ind.x_va_ind)';
    x_out.x_qa = x(:,Sim.Ind.x_qa_ind)';
    x_out.x_pa = x(:,Sim.Ind.x_pa_ind)';
end