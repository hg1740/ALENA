function [Matrices,Aero,MassProps] = InitialiseRobbie_old(aircraft_data,extra_flag,hp)

if nargin < 2 || isempty(extra_flag)
    extra_flag = 0;
end
if nargin < 3
    hp         = [];
end

%REMOVE THIS HEINOUS USE OF MATLAB CODE FOREVER!
%   - Robbie is calling a script which then populates workspace variables
%   that are then used in subsequent analysis. This is a poor (yet
%   practical) use of MATLAB. Better to make 'ReadDario' a function!
% ReadDario(aircraft_data)
ReadDario

if (extra_flag == 1) && (n_elem > 1)
    fprintf('\nThis method is not supported for aircraft with more than one member. Only using the first defined member for the analysis...\n')
    n_elem = n_elem(1);
end

Dtot    = cell(length(n_elem),1); %call length(n_elem) once and save as a variable. Will be slightly quicker.
Dtot2   = cell(length(n_elem),1);
Btot    = cell(length(n_elem),1);
Btot2   = cell(length(n_elem),1);
Btot_mp = cell(length(n_elem),1);
Btot_mp_pnt = cell(length(n_elem),1);
B_aero  = cell(length(n_elem),1);
e1tot   = cell(length(n_elem),1);
e1tot2  = cell(length(n_elem),1);
Mtot    = cell(length(n_elem),1);
Mtot2   = cell(length(n_elem),1);
Ftot    = cell(length(n_elem),1);
Ftot2   = cell(length(n_elem),1);
Gvtot   = cell(length(n_elem),1);
Gftot   = cell(length(n_elem),1);
G_bc    = cell(length(n_elem),1);

R1_bc   = cell(length(n_elem),1);
R2_bc   = cell(length(n_elem),1);
R3_bc   = cell(length(n_elem),1);
R4_bc   = cell(length(n_elem),1);

loadInt     = cell(length(n_elem),1);
loadInt2    = cell(length(n_elem),1);
loadIntAero = cell(length(n_elem),1);
PntIn       = cell(length(n_elem),1);

Cfull         = cell(length(n_elem),1);
CFFbar_full   = cell(length(n_elem),1);

Ab_sep  = cell(length(n_elem),1);
Ab2_sep = cell(length(n_elem),1);
Aa_sep  = cell(length(n_elem),1);
Aa2_sep = cell(length(n_elem),1);

Ab1var_sep = cell(length(n_elem),1);
Ab2var_sep = cell(length(n_elem),1);
Ab3var_sep = cell(length(n_elem),1);
Ab4var_sep = cell(length(n_elem),1);
Aa1var_sep = cell(length(n_elem),1);
Aa2var_sep = cell(length(n_elem),1);
Aa3var_sep = cell(length(n_elem),1);
Aa4var_sep = cell(length(n_elem),1);
Aa5var_sep = cell(length(n_elem),1);
Aa6var_sep = cell(length(n_elem),1);

e1tilde = skew([1;0;0]);
e2tilde = skew([0;1;0]);
e3tilde = skew([0;0;1]);
E1tilde = [ zeros(3), -e1tilde,   zeros(3),-e2tilde,  zeros(3), -e3tilde,   zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3);
           -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3),  zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde];
E2tilde = [ zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3), -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3);
           -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3),  zeros(3), -e1tilde,   zeros(3),-e2tilde,  zeros(3), -e3tilde];
E3tilde = [ zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde,  -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3);
            zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3),  zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde];
P_36    = sparse(36,36);
for ii = 1:6
    ind6       = [1:6] + (ii-1)*6;
    for kk = 1:6
        ind6_2 = [1:6] + (kk-1)*6;
        P_36(ind6(kk),ind6_2(ii)) = 1;
    end
end
        
% if ~exist('p0','var') || isempty(p0{jj})
%     p0     = cell(1,length(n_elem));
% end
if ~exist('eps0','var') %TODO - Try to avoid using 'exist' where possible.
    eps0   = cell(1,length(n_elem));
end
if ~exist('kappa0','var')
    kappa0 = cell(1,length(n_elem));
end
if ~exist('CBB','var')
    CBB = cell(1,length(n_elem));
end
if ~exist('Parent','var')
    Parent = cell(1,length(n_elem));
end
if ~exist('CFFbar','var')
    CFFbar = cell(1,length(n_elem));
end
if ~exist('M_pt_store','var')
    M_pt_store = cell(1,length(n_elem));
end

Children = cell(1,length(n_elem));

for jj = 1:length(n_elem)
    Cfull_elem  = [];
    CFFbar_elem = [];     
    
    Dtot{jj}    = zeros((n_elem(jj)+1)*6,n_elem(jj)*6);
    Dtot2{jj}   = zeros( n_elem(jj)*12  ,n_elem(jj)*12);
    Btot{jj}    = zeros((n_elem(jj)+1)*6);
    Btot2{jj}   = zeros( n_elem(jj)*12  ,n_elem(jj)*12);
    Btot_mp{jj} = zeros((n_elem(jj)+1)*6,n_elem(jj)*6);
    Btot_mp_pnt{jj} = zeros((n_elem(jj)+1)*6,n_elem(jj)*6);
    B_aero{jj}  = zeros((n_elem(jj)+1)*6,length(s_out{jj})*6);
    e1tot{jj}   = zeros( n_elem(jj)*6   ,1);
    e1tot2{jj}  = zeros( n_elem(jj)*12  ,1);
    Gvtot{jj}    = zeros((n_elem(jj)+1)*12);
    Gftot{jj}    = zeros((n_elem(jj)+1)*12);
    G_bc{jj}    = zeros( n_elem(jj)*12  ,6);
    
    G_bc{jj}(1:6,:) = eye(6);
    
    if strcmp(type{jj},'LS')
        isLS{jj} = 1;
    elseif strcmp(type{jj},'FUS')
        isLS{jj} = 0;
    else
        isLS{jj} = 0;
    end

    R1_bc{jj}            = sparse(zeros(n_elem(jj)*6,(n_elem(jj)+1)*6));
    R1_bc{jj}(:,7:end)   = eye(n_elem(jj)*6);
    R2_bc{jj}            = sparse(zeros(n_elem(jj)*6,(n_elem(jj)+1)*6));
    R2_bc{jj}(:,1:end-6) = eye(n_elem(jj)*6);

    if ~exist('p0','var')     || isempty(p0{jj})
        p0{jj} = [0;0;0];
    end
    if ~exist('eps0','var')   || isempty(eps0{jj})
        eps0{jj}   = zeros(3,1,n_elem(jj));
    end
    if ~exist('kappa0','var') || isempty(kappa0{jj})
        kappa0{jj} = zeros(3,1,n_elem(jj));
    end
    if ~exist('CBB','var')    || isempty(CBB{jj})
        CBB{jj} = repmat(eye(3),1,1,n_elem(jj));
    end
    if ~exist('Parent','var') || isempty(Parent{jj})
        Parent{jj} = 0;
    end
    if ~exist('CFFbar','var') || isempty(CFFbar{jj})
        CFFbar{jj} = repmat(eye(6),1,1,n_elem(jj));
    end
    if isempty(M_pt_store{jj})
        M_pt_store{jj} = zeros(6,6,n_elem(jj)+1);
    end
    
    if Parent{jj}(1)
        Children{Parent{jj}(1)} = [Children{Parent{jj}(1)};jj];
        parent_ind               = Parent{jj}(1);
        ind_tmp                  = [1:6] + (Parent{jj}(2)-2)*6;
        R3_bc{jj}                = sparse(zeros(n_elem(parent_ind)*6,(n_elem(jj)+1)*6));
        R3_bc{jj}(ind_tmp,1:6)   = blkdiag(CaB0{Parent{jj}(1)}'*CaB0{jj},CaB0{Parent{jj}(1)}'*CaB0{jj});
        R4_bc{jj}                = sparse(zeros(n_elem(jj)*6,(n_elem(parent_ind)+1)*6));
        R4_bc{jj}(1:6,end-5:end) = blkdiag(CaB0{Parent{jj}(1)}'*CaB0{jj},CaB0{Parent{jj}(1)}'*CaB0{jj})'; 
    end
    
    vel_in{jj}          = zeros((n_elem(jj)+1)*6,6);
    vel_in{jj}(1:6,1:6) = [CaB0{jj}',zeros(3);zeros(3),CaB0{jj}']*[eye(3),-skew(p0{jj});zeros(3),eye(3)];
    
    Mtot{jj}    = zeros((n_elem(jj)+1)*6, n_elem(jj)*6 );
    Mtot2{jj}   = zeros( n_elem(jj)   *12,n_elem(jj)*12);
    Ftot{jj}    = zeros((n_elem(jj)+1)*6, n_elem(jj)*6 );
    Ftot2{jj}   = zeros( n_elem(jj)   *12,n_elem(jj)*12);
    
    loadInt{jj}     = zeros(6,(n_elem(jj)+1)*6);
    loadInt2{jj}    = zeros(6, n_elem(jj)   *12);
    loadIntAero{jj} = zeros(6,length(s_out{jj})*6);
    PntInt{jj}      = zeros(6, n_elem(jj)   *12);
    
    B1 = zeros(n_elem(jj)*3+6,n_elem(jj)*3+3);
    B1(1:end-3,:) = eye(n_elem(jj)*3+3);
    B1(4:end,:)   = B1(4:end,:) + eye(n_elem(jj)*3+3);
    B1 = B1(1:end-3,:);
    
    B2 = zeros(n_elem(jj)*3+3,n_elem(jj)*3);
    B2(4:end,:)   = 2*eye(n_elem(jj)*3);
    
    x_p_recov1{jj} = inv(B1)*B2;
    x_p_recov2{jj} = zeros(n_elem(jj)*3+3,3);
    x_p_recov2{jj}(1:3,:) = eye(3);
    x_p_recov2{jj} = inv(B1)*x_p_recov2{jj};
    
    B1 = zeros(n_elem(jj)*6+12,n_elem(jj)*6+6);
    B1(1:end-6,:) = eye(n_elem(jj)*6+6);
    B1(7:end,:)   = B1(7:end,:) + eye(n_elem(jj)*6+6);
    B1 = B1(1:end-6,:);
    
    B2 = zeros(n_elem(jj)*6+6,n_elem(jj)*6);
    B2(7:end,:)   = 2*eye(n_elem(jj)*6);
    
    x_v_recov1_tmp = inv(B1)*B2;
    x_v_recov2_tmp = zeros(n_elem(jj)*6+6,6);
    x_v_recov2_tmp(1:6,:) = eye(6);
    x_v_recov2_tmp = inv(B1)*x_v_recov2_tmp;
    
    x_v_recov1{jj} = [];
    x_v_recov2{jj} = [];
    for ii = 1:n_elem(jj)+1
        ind = [1:6] + (ii-1)*6;
        orig_v         = zeros(6,size(x_v_recov1_tmp,2));
        orig_v(:,ind)  = eye(6);
        if ii == n_elem(jj)+1
            x_v_recov1{jj} = [x_v_recov1{jj};x_v_recov1_tmp(ind,:)];
            x_v_recov2{jj} = [x_v_recov2{jj}();x_v_recov2_tmp(ind,:)];
        else
            x_v_recov1{jj} = [x_v_recov1{jj};x_v_recov1_tmp(ind,:);orig_v];
            x_v_recov2{jj} = [x_v_recov2{jj}();x_v_recov2_tmp(ind,:);zeros(6)];
        end
    end
    
    if extra_flag
        Q1_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q2_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q3_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q4_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q5_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q6_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q7_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q8_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)*6    ));
        Q9_quad{jj}  = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q10_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q11_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q12_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q13_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q14_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q15_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        Q16_quad{jj} = sparse(zeros((n_elem(jj)+1)*6,  n_elem(jj)^2*36 ));
        P_quad       = sparse(       n_elem(jj)^2*144, n_elem(jj)^2*144);
%         P_quad_test  = sparse(       n_elem(jj)^2*144, n_elem(jj)^2*144);
    end
    
    ind5_start  = 1;
    for ii = 1:n_elem(jj)       
        M = [m{jj}(ii)*eye(3)              ,-m{jj}(ii)*skew(cg{jj}(:,:,ii));
             m{jj}(ii)*skew(cg{jj}(:,:,ii)), eye(3)*J{jj}(:,:,ii)   ];
        
        Cfull_elem  = blkdiag(Cfull_elem, C{jj}(:,:,ii));
        CFFbar_elem = blkdiag(CFFbar_elem,CFFbar{jj}(:,:,ii));
        
        e1 = zeros(6,1); e1(1) = 1;
        e1 = e1 + [eps0{jj}(:,:,ii);kappa0{jj}(:,:,ii)];
        e1 = K{jj}(:,:,ii)*e1;
        
        ind1 = [1:12] + (ii-1)*12;
        ind2 = [1:12] + (ii-1)*6;
        ind3 = [1:6]  + (ii-1)*6;
                
        CBB_tot = blkdiag(CBB{jj}(:,:,ii),CBB{jj}(:,:,ii),eye(6));
        
        A1 = [-eye(6)*(s{jj}(ii+1)-s{jj}(ii)),  zeros(6)]/(s{jj}(ii)-s{jj}(ii+1))*CBB_tot;
        A2 = [ eye(6)                        , -eye(6)  ]/(s{jj}(ii)-s{jj}(ii+1))*CBB_tot;
        
        Ab  = A1'*(s{jj}(ii+1)-s{jj}(ii)) + A2'*(s{jj}(ii+1)-s{jj}(ii))^2/2;
        Ab2 = eye(6)*(s{jj}(ii+1)-s{jj}(ii));
        Aa  = C{jj}(:,:,ii);
        Aa2 = M;
        
        %%
        Ab1var  = (A1'*A1*(s{jj}(ii+1)-s{jj}(ii))     + A1'*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2 + ...
                   A2'*A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A2'*A2*(s{jj}(ii+1)-s{jj}(ii))^3/3);
        
        Ab2var  = (A1'*A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A1'*A2*(s{jj}(ii+1)-s{jj}(ii))^3/3 + ...
                   A2'*A1*(s{jj}(ii+1)-s{jj}(ii))^3/3 + A2'*A2*(s{jj}(ii+1)-s{jj}(ii))^4/4);
        
        Ab3var  = (A1*(s{jj}(ii+1)-s{jj}(ii))     + A2*(s{jj}(ii+1)-s{jj}(ii))^2/2);
        
        Ab4var  = (A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A2*(s{jj}(ii+1)-s{jj}(ii))^3/3);
        
        Aa1var  = C{jj}(:,:,ii)*A1;
        Aa2var  = C{jj}(:,:,ii)*A2;
        Aa3var  =             M*A1;
        Aa4var  =             M*A2;
        Aa5var  =               A1;
        Aa6var  =               A2;
        %%
        
        Bdst = A1'*A1*(s{jj}(ii+1)-s{jj}(ii))     + A1'*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2 + ...
               A2'*A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A2'*A2*(s{jj}(ii+1)-s{jj}(ii))^3/3;
        
        D   = -A2'*eye(6)*(s{jj}(ii+1)-s{jj}(ii));
%         D2  = -(A2'*A1*(s{jj}(ii+1)-s{jj}(ii)) + A2'*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2);
        D2  = (A1'*A2*(s{jj}(ii+1)-s{jj}(ii)) + A2'*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2);
        
        Ms  = Ab  * Aa2;
        Ms2 = A1'*M*A1*(s{jj}(ii+1)-s{jj}(ii))     + A1'*M*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2 + ...
              A2'*M*A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A2'*M*A2*(s{jj}(ii+1)-s{jj}(ii))^3/3;
        Fs  = Ab  * Aa;
        Fs2 = A1'*C{jj}(:,:,ii)*A1*(s{jj}(ii+1)-s{jj}(ii))     + A1'*C{jj}(:,:,ii)*A2*(s{jj}(ii+1)-s{jj}(ii))^2/2 + ...
              A2'*C{jj}(:,:,ii)*A1*(s{jj}(ii+1)-s{jj}(ii))^2/2 + A2'*C{jj}(:,:,ii)*A2*(s{jj}(ii+1)-s{jj}(ii))^3/3;
        
        C_dihedral_v = [cos(local_dihedral{jj}(ii)) 0 -sin(local_dihedral{jj}(ii)) ; 0 1 0 ; sin(local_dihedral{jj}(ii)) 0 cos(local_dihedral{jj}(ii))];
        C_sweep_v    = [cos(local_sweep{jj}(ii)) -sin(local_sweep{jj}(ii)) 0 ; sin(local_sweep{jj}(ii)) cos(local_sweep{jj}(ii)) 0 ; 0 0 1];
        C_twist_v    = [1 0 0 ; 0 cos(local_twist{jj}(ii)) sin(local_twist{jj}(ii)) ; 0 -sin(local_twist{jj}(ii)) cos(local_twist{jj}(ii))];
          
        C_rot_v      = C_dihedral_v*C_sweep_v*C_twist_v;
        
        C_dihedral_f = [cos(local_dihedral{jj}(ii+1)) 0 -sin(local_dihedral{jj}(ii+1)) ; 0 1 0 ; sin(local_dihedral{jj}(ii+1)) 0 cos(local_dihedral{jj}(ii+1))];
        C_sweep_f    = [cos(local_sweep{jj}(ii+1)) -sin(local_sweep{jj}(ii+1)) 0 ; sin(local_sweep{jj}(ii+1)) cos(local_sweep{jj}(ii+1)) 0 ; 0 0 1];
        C_twist_f    = [1 0 0 ; 0 cos(local_twist{jj}(ii+1)) sin(local_twist{jj}(ii+1)) ; 0 -sin(local_twist{jj}(ii+1)) cos(local_twist{jj}(ii+1))];
          
        C_rot_f      = C_dihedral_f*C_sweep_f*C_twist_f;
        
        C_rot{jj}(:,:,ii) = C_dihedral_v*C_sweep_v*C_twist_v';
        if ii == n_elem(jj)
            C_rot{jj}(:,:,ii+1) = C_dihedral_f*C_sweep_f*C_twist_f';
        end
        
%         G   = [zeros(6),zeros(6);eye(6),-eye(6)];
        Gv   = [zeros(6),zeros(6);eye(6),-blkdiag(C_rot_v,C_rot_v)];
        Gf   = [zeros(6),zeros(6);eye(6),-blkdiag(C_rot_f,C_rot_f)];
        
        nodes_in_elem = find(s_out{jj}>=s{jj}(ii)&s_out{jj}<=s{jj}(ii+1));
        B_aero_e      = zeros(12,length(nodes_in_elem)*6);
        loadIntAero_e = zeros( 6,length(nodes_in_elem)*6);
        for kk = 1:length(nodes_in_elem)-1
            ind4 = [1:12] + (kk-1)*6;
            ds_elem    = s_out{jj}(nodes_in_elem(kk+1))-s_out{jj}(nodes_in_elem(kk));
            B_aero_e_n = A1'*A1*(ds_elem)     + A1'*A2*(ds_elem)^2/2 + ...
                         A2'*A1*(ds_elem)^2/2 + A2'*A2*(ds_elem)^3/3;
                     
            B_aero_e(:,ind4)      = B_aero_e(:,ind4) + B_aero_e_n;
            loadIntAero_e(:,ind4) = loadIntAero_e(:,ind4) + Ab';
        end
        ind5       = ind5_start:ind5_start+length(nodes_in_elem)*6-1;
        ind5_start = ind5_start+length(nodes_in_elem)*6-6;
        
        %% Alt cross method (using kron)
        MkronI                     = kron(M,                        eye(6));
        CkronI                     = kron(C{jj}(:,:,ii),            eye(6));
        Ikronef                    = kron(eye(6),                   e1    );
        MptkronI                   = kron(M_pt_store{jj}(:,:,ii+1), eye(6));
        
        e1ds                       = [1;0;0]*(s{jj}(ii+1)-s{jj}(ii))/2;
        T_pt{jj}(:,:,ii)           = [eye(3) -skew(e1ds);zeros(3) eye(3)];
        T_ptkronT_pt               = kron(T_pt{jj}(:,:,ii),T_pt{jj}(:,:,ii));
        
        C_pre1{jj}(:,:,ii)         = Ab*E1tilde*MkronI;
        C_pre2{jj}(:,:,ii)         = C_pre1{jj}(:,:,ii)*P_36;
        A_pre1{jj}(:,:,ii)         = Ab*E2tilde*CkronI;
        A_pre2{jj}(:,:,ii)         = A_pre1{jj}(:,:,ii)*P_36;
        A_pre3{jj}(:,:,ii)         = A_pre2{jj}(:,:,ii)*Ikronef;
        E_pre1{jj}(:,:,ii)         = Ab*E3tilde*CkronI;
        E_pre2{jj}(:,:,ii)         = E_pre1{jj}(:,:,ii)*P_36;
        E_pre3{jj}(:,:,ii)         = E_pre2{jj}(:,:,ii)*Ikronef;
        E2tildeIkronM{jj}(:,:,ii)  = E2tilde*kron(eye(6),M);
        
        C_pt_pre1{jj}(:,:,ii)      = E1tilde*MptkronI*T_ptkronT_pt;
        C_pt_pre2{jj}(:,:,ii)      = C_pt_pre1{jj}(:,:,ii)*P_36;
        
        %% Combine Everything
        Ab_sep{jj}(:,:,ii)  = Ab;
        Ab2_sep{jj}(:,:,ii) = Ab2;
        Aa_sep{jj}(:,:,ii)  = Aa;
        Aa2_sep{jj}(:,:,ii) = Aa2;
        
        Ab1var_sep{jj}(:,:,ii)  = Ab1var;
        Ab2var_sep{jj}(:,:,ii)  = Ab2var;
        Ab3var_sep{jj}(:,:,ii)  = Ab3var;
        Ab4var_sep{jj}(:,:,ii)  = Ab4var;
        Aa1var_sep{jj}(:,:,ii)  = Aa1var;
        Aa2var_sep{jj}(:,:,ii)  = Aa2var;
        Aa3var_sep{jj}(:,:,ii)  = Aa3var;
        Aa4var_sep{jj}(:,:,ii)  = Aa4var;
        Aa5var_sep{jj}(:,:,ii)  = Aa5var;
        Aa6var_sep{jj}(:,:,ii)  = Aa6var;
        
        Dtot{jj}(ind2,ind3)        = Dtot{jj}(ind2,ind3)     + D;
        Dtot2{jj}(ind1,ind1)       = Dtot2{jj}(ind1,ind1)    + D2;
        Btot{jj}(ind2,ind2)        = Btot{jj}(ind2,ind2)     + Bdst;%
        Btot2{jj}(ind1,ind1)       = Btot2{jj}(ind1,ind1)    + Bdst;%
        Btot_mp{jj}(ind2,ind3)     = Btot_mp{jj}(ind2,ind3)  + Ab;%
        Btot_mp_pnt{jj}(ind2,ind3) = Btot_mp_pnt{jj}(ind2,ind3)  + [eye(6);eye(6)]/2;%
        B_aero{jj}(ind2,ind5)      = B_aero{jj}(ind2,ind5)   + B_aero_e;
        e1tot{jj}(ind3,1)          = e1;
        e1tot2{jj}(ind1,1)         = [e1;e1];
        Mtot{jj}(ind2,ind3)        = Mtot{jj}(ind2,ind3)      + Ms;
        Mtot2{jj}(ind1,ind1)       = Mtot2{jj}(ind1,ind1)     + Ms2;% + M_pnt{jj}(ind1,ind1);
        Ftot{jj}(ind2,ind3)        = Ftot{jj}(ind2,ind3)      + Fs;
        Ftot2{jj}(ind1,ind1)       = Ftot2{jj}(ind1,ind1)     + Fs2;
        Gvtot{jj}(ind1+6,ind1+6)   = Gvtot{jj}(ind1+6,ind1+6) + Gv;
        Gftot{jj}(ind1+6,ind1+6)   = Gftot{jj}(ind1+6,ind1+6) + Gf;
        
        Mtot{jj}(ind3+6,ind3)      = Mtot{jj}(ind3+6,ind3)     + M_pt_store{jj}(:,:,ii+1)*T_pt{jj}(:,:,ii);
        
        loadInt{jj}(1:6,ind2)      = loadInt{jj}(1:6,ind2)     + Ab';
        loadInt2{jj}(1:6,ind1)     = loadInt2{jj}(1:6,ind1)    + Ab';
        loadIntAero{jj}(1:6,ind5)  = loadIntAero{jj}(1:6,ind5) + loadIntAero_e;
        PntInt{jj}(1:6,ind1)       = [eye(6),eye(6)];

        if extra_flag
            ind36 = [1:36] + (ii-1)*36;
            Q1_quad{jj}( ind2,ind3 ) = Q1_quad{jj}( ind2,ind3 ) + Ms;
            Q2_quad{jj}( ind2,ind3 ) = Q2_quad{jj}( ind2,ind3 ) + D*1e-6;
            Q4_quad{jj}( ind2,ind3 ) = Q4_quad{jj}( ind2,ind3 ) + Fs;
            Q6_quad{jj}( ind2,ind3 ) = Q6_quad{jj}( ind2,ind3 ) - D;
            Q7_quad{jj}( ind2,ind3 ) = Q7_quad{jj}( ind2,ind3 ) - D;
            Q9_quad{jj}( ind2,ind36) = Q9_quad{jj}( ind2,ind36) + Ab*E1tilde*kron(M,            eye(6));
            Q12_quad{jj}(ind2,ind36) = Q12_quad{jj}(ind2,ind36) + Ab*E2tilde*kron(C{jj}(:,:,ii),eye(6));
            Q14_quad{jj}(ind2,ind36) = Q14_quad{jj}(ind2,ind36) + Ab*E3tilde*kron(C{jj}(:,:,ii),eye(6));
        end
    end
    if extra_flag
        for ii = 1:(12*n_elem(jj))
            ind192       = [1:12*n_elem(jj)] + (ii-1)*12*n_elem(jj);
            for kk = 1:(12*n_elem(jj))
                ind192_2 = [1:12*n_elem(jj)] + (kk-1)*12*n_elem(jj);
                P_quad(ind192(kk),ind192_2(ii)) = 1;
            end
        end
    end
    Cfull{jj}       = Cfull_elem;
    CFFbar_full{jj} = CFFbar_elem;
end

if ~exist('CBA0','var')
    CBA0 = [];
end
if ~exist('le','var')
    le = [];
end
if ~exist('te','var')
    te = [];
end
if ~exist('cp','var')
    cp = [];
end
if ~exist('RenderPoints','var')
    RenderPoints = [];
end
if ~exist('RenderPoints_RB','var')
    RenderPoints_RB = [];
end
if ~exist('dCLdd','var')
    dCLdd = [];
end
if ~exist('dCMdd','var')
    dCMdd = [];
end
if ~exist('delta_flap','var')
    delta_flap = [];
end

Matrices.Ab  = Ab_sep;
Matrices.Ab2 = Ab2_sep;
Matrices.Aa  = Aa_sep;
Matrices.Aa2 = Aa2_sep;

Matrices.Ab1var = Ab1var_sep;
Matrices.Ab2var = Ab2var_sep;
Matrices.Ab3var = Ab3var_sep;
Matrices.Ab4var = Ab4var_sep;
Matrices.Aa1var = Aa1var_sep;
Matrices.Aa2var = Aa2var_sep;
Matrices.Aa3var = Aa3var_sep;
Matrices.Aa4var = Aa4var_sep;
Matrices.Aa5var = Aa5var_sep;
Matrices.Aa6var = Aa6var_sep;

Matrices.C_pre1    = C_pre1;
Matrices.C_pre2    = C_pre2;
Matrices.A_pre1    = A_pre1;
Matrices.A_pre2    = A_pre2;
Matrices.A_pre3    = A_pre3;
Matrices.E_pre1    = E_pre1;
Matrices.E_pre2    = E_pre2;
Matrices.E_pre3    = E_pre3;
Matrices.C_pt_pre1 = C_pt_pre1;
Matrices.C_pt_pre2 = C_pt_pre2;
Matrices.E1tilde = E1tilde;
% Matrices.E2tilde = E2tilde;
Matrices.E2tildeIkronM  = E2tildeIkronM;

Matrices.Dtot    = Dtot;
Matrices.Dtot2   = Dtot2;
Matrices.Btot    = Btot;
Matrices.Btot2   = Btot2;
Matrices.Btot_mp = Btot_mp;
Matrices.Btot_mp_pnt = Btot_mp_pnt;
Matrices.B_aero  = B_aero;

Matrices.Gvtot    = Gvtot;
Matrices.Gftot    = Gftot;
Matrices.G_bc     = G_bc;

Matrices.R1_bc    = R1_bc;
Matrices.R2_bc    = R2_bc;
Matrices.R3_bc    = R3_bc;
Matrices.R4_bc    = R4_bc;

Matrices.vel_in   = vel_in;

Matrices.CFFbar   = CFFbar_full;

Matrices.M   = Mtot;
Matrices.M2  = Mtot2;%M_pnt;%
Matrices.T1  = Ftot;
Matrices.T12 = Ftot2;

Matrices.e1tot  = e1tot;
Matrices.e1tot2 = e1tot2;

Matrices.eps0   = eps0;
Matrices.kappa0 = kappa0;

Matrices.M_pt = M_pt_store;
Matrices.T_pt = T_pt;

Matrices.m  = m;
Matrices.J  = J;
Matrices.cg = cg;

Matrices.C = C;

Matrices.C_rot = C_rot;

Matrices.loadInt     = loadInt;
Matrices.loadInt2    = loadInt2;
Matrices.loadIntAero = loadIntAero;
Matrices.pntInt      = PntInt;

Matrices.Cfull = Cfull;

Matrices.n_elem = n_elem;
Matrices.l      = l;

Matrices.type   = type;

Matrices.s             = s;
Matrices.s2            = s2;
Matrices.s_node_ind    = s_node_ind;
Matrices.s_node_ind2   = s_node_ind2;
Matrices.s_aero        = s_aero;
Matrices.s_aero_mp     = s_aero_mp;
Matrices.s_out         = s_out;
Matrices.s_out2        = s_out2;
Matrices.s_mp_ind      = s_mp_ind;
Matrices.s_aero_ind    = s_aero_ind;
Matrices.s_aero_mp_ind = s_aero_mp_ind;

Matrices.p0            = p0;
Matrices.CaB0          = CaB0;
Matrices.CBA0          = CBA0;
Matrices.CBB           = CBB;

Matrices.chord         = chord;
Matrices.le            = le;
Matrices.te            = te;
Matrices.cp            = cp;
% Matrices.xs            = xs;

Matrices.RenderPoints     = RenderPoints;
Matrices.RenderPoints_RB  = RenderPoints_RB;

Matrices.x_p_recov1    = x_p_recov1;
Matrices.x_p_recov2    = x_p_recov2;
Matrices.x_v_recov1    = x_v_recov1;
Matrices.x_v_recov2    = x_v_recov2;

Matrices.Parent        = Parent;
Matrices.Children      = Children;

if extra_flag
    reorder_mat1            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat2            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat3            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat4            = sparse(36*n_elem^2,144*n_elem^2);
    
    counter = 0;
    for ii = 1:n_elem
        ind6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            VO_tmp1(ind6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                VO_tmp2(ind6(kk))       = 1;
                reorder_mat1(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    counter = 0;
    for ii = 1:n_elem
        ind6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1(ind6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                VO_tmp2(ind6(kk))       = 1;
                reorder_mat2(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    counter = 0;
    for ii = 1:n_elem
        ind6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1(ind6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2(ind6(kk))       = 1;
                reorder_mat4(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    reorder_mat             = [reorder_mat1;reorder_mat2;reorder_mat3;reorder_mat4];
    reorder_mat_red         = reorder_mat;
    reorder_mat_red(sum(reorder_mat,2)==0,:) = [];
    
    Q1_tot                  = [R1_bc{1}*Q1_quad{1}, R1_bc{1}*Q2_quad{1} ;...
                               R2_bc{1}*Q3_quad{1}, R2_bc{1}*Q4_quad{1} ];
    Q2_tot                  = [R1_bc{1}*Q5_quad{1}, R1_bc{1}*Q6_quad{1} ;...
                               R2_bc{1}*Q7_quad{1}, R2_bc{1}*Q8_quad{1} ];
    Q3_tot                  = [R1_bc{1}*Q9_quad{1}, R1_bc{1}*Q10_quad{1},R1_bc{1}*Q11_quad{1},R1_bc{1}*Q12_quad{1};...
                               R2_bc{1}*Q13_quad{1},R2_bc{1}*Q14_quad{1},R2_bc{1}*Q15_quad{1},R2_bc{1}*Q16_quad{1}];
    Matrices.A_quad         = Q2_tot;%Q1_tot\Q2_tot;
    Matrices.B_quad         = Q3_tot*reorder_mat;%Q1_tot\Q3_tot*reorder_mat;
    Matrices.C_quad         = Q1_tot;
    Matrices.I_quad         = speye(size(Q1_tot));
    Matrices.P_quad         = P_quad;
    Matrices.B_x_P_quad_p_I = Matrices.B_quad*(P_quad + speye(size(P_quad)));
    Matrices.e1_quad        = sparse([zeros(size(e1tot{1}));e1tot{1}]);
    Matrices.A_quad         = Matrices.A_quad + Matrices.B_quad*kron(Matrices.e1_quad,Matrices.I_quad);
    Matrices.D_quad         = sparse(n_elem*12,(n_elem+1)*6);
    Matrices.D_quad(1:(n_elem+1)*6,1:(n_elem+1)*6) = speye((n_elem+1)*6);
end

%%
Aero.chord = chord;
Aero.xs    = xs;
Aero.tau   = tau;

Aero.lift_dist = lift_dist;

Aero.dCLdd = dCLdd;
Aero.dCMdd = dCMdd;
Aero.delta_flap = delta_flap;

Matrices.isLS  = isLS;

%% Initialise labels
loadLabels  = {'\epsilon_x','\epsilon_y','\epsilon_z','\kappa_x (m^{-1})','\kappa_y (m^{-1})','\kappa_z (m^{-1})'};
forceLabels = {'F_x (N)','F_y (N)','F_z (N)','M_x (Nm)','M_y (Nm)','M_z (Nm)'};

%% Calculate Mass Properties
MassProps = calculateCoM(x,Matrices,1);

%% Change Mass
% deltaMass              = 61441.5892 - MassProps.m;
% Matrices.M_rb(1:3,1:3) = Matrices.M_rb(1:3,1:3) + eye(3)*deltaMass;
% MassProps = calculateCoM(x,Matrices,1);

%% Clear all variables from the Workspace, but for x, Matrices, Aero and MassProps
% clearvars -except x Matrices MassProps Aero ff flex_vec