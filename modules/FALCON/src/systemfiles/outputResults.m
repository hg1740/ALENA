function OUT = outputResults(X,SYS,Matrices,r_full,Sim)

if Sim.Soln == 1
    OUT.total_Aero   = zeros(6,1);
    OUT.total_Grav   = zeros(6,1);
    OUT.total_Aero_G = zeros(6,1);
    OUT.total_Grav_G = zeros(6,1);
    OUT.tail_Aero_G  = zeros(6,1);
    OUT.wing_Aero_G  = zeros(6,1);
    CGa          = Matrices.CGa;
    for jj = 1:length(Matrices.n_elem)
        OUT.x_f{jj}    = X.x_f_p1{jj};
        [~ ,CaB ] = strains2coords_all(X.x_f_p1,Matrices,jj);
        
        OUT.aeroForce{jj}   = SYS.f{jj};
        OUT.Local_AoA{jj}     = SYS.AoA_qs_t{jj}*180/pi;
        
        for ii = 1:Matrices.n_elem(jj)
            ind3  = [1:6]  + (ii-1)*6;
            CGB_e = Matrices.CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii));
            OUT.aeroForce_G{jj}(ind3,1) = [CGB_e,zeros(3);zeros(3),CGB_e]*SYS.f{jj}(ind3);
        end
        
        OUT.total_Aero   = OUT.total_Aero   + SYS.Btot_r{jj}*SYS.f{jj};
        OUT.total_Grav   = OUT.total_Grav   + SYS.Btot_r{jj}*SYS.f_grav{jj};
        OUT.total_Aero_G = OUT.total_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
        OUT.total_Grav_G = OUT.total_Grav_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f_grav{jj};
        if jj == 1 || jj == 2
            OUT.wing_Aero_G = OUT.wing_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
        elseif jj == 3 || jj == 4
            OUT.tail_Aero_G = OUT.tail_Aero_G + [CGa zeros(3);zeros(3) CGa]*SYS.Btot_r{jj}*SYS.f{jj};
        end
        OUT.res = norm(r_full);
    end
    OUT.total_Grav   = OUT.total_Grav   +                              Matrices.M_rb*[CGa' zeros(3); zeros(3) CGa']*Sim.grav_vec;
    OUT.total_Grav_G = OUT.total_Grav_G + [CGa zeros(3); zeros(3) CGa]*Matrices.M_rb*[CGa' zeros(3); zeros(3) CGa']*Sim.grav_vec;
elseif Sim.Soln == 2
    OUT.x_f     = X.x_f;
    OUT.x_v     = X.x_v;
    OUT.x_q     = X.x_q;
    OUT.x_p     = X.x_p;
    OUT.x_x     = X.x_x;
    
    OUT.x_f_dot = X.x_f_dot;
    OUT.x_v_dot = X.x_v_dot;
    OUT.x_q_dot = X.x_q_dot;
    OUT.x_p_dot = X.x_p_dot;
    OUT.x_x_dot = X.x_x_dot;
    
    OUT.x_va = X.x_va;
    OUT.x_qa = X.x_qa;
    OUT.x_pa = X.x_pa;
    
    OUT.x_vg = X.x_vg;
    
    OUT.x_va_dot = X.x_va_dot;
    OUT.x_qa_dot = X.x_qa_dot;
    OUT.x_pa_dot = X.x_pa_dot;
    
%     OUT.KE      = X.KE;
%     OUT.PE      = X.PE;
    
    OUT.res     = X.residual;
    
    OUT.aeroForce     = X.Force_Store;
    OUT.total_Aero    = X.total_Aero;
    OUT.total_Aero_G  = X.total_Aero_G;
    OUT.wing_Aero_G   = X.wing_Aero_G;
    OUT.tail_Aero_G   = X.tail_Aero_G;
    OUT.Local_AoA     = X.Local_AoA;
    OUT.RBEulerAngles = X.RBEulerAngles;
end