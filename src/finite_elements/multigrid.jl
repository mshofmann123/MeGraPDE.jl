# © Anna Weller, University of Cologne, 2023

using DocStringExtensions, SparseArrays

"""
    SmoothingMatrices for weighted Jacobi Iterations

$(FIELDS)
"""
struct SmoothingMatrices
    "D^{-1} where D is diagonal"
    D_inv:: Matrix
    "Lower Triangle Matrix"
    Low:: Matrix
    "Upper Triangle Matrix"
    Up:: Matrix
end

"""
    smooth(u::Vector, b::Vector, smooth_mats::SmoothingMatrices)

Perform one smoothing step of the weighted jacobi iteration u^(k+1) = 1/2*(I-D^{-1}*(Low+Up))*u + 1/2*D^{-1}*b. 

Objective is the iterative solution of Ax=b.

# Arguments 
- 'u::Vector': current iterative  
- 'b::Vector': right-hand side
- 'smooth_mats::SmoothingMatrices': Matrices smooth_mats.D_inv,smooth_mats.Low,smooth_mats.Up of the decomposition A=Low+Up+D 

"""
function smooth(u::Vector, b::Vector, smooth_mats::SmoothingMatrices)
    return 0.5*(I-smooth_mats.D_inv*(smooth_mats.Low+smooth_mats.Up))*u+0.5*smooth_mats.D_inv*(b);
end



function P_J(u_tp1,u_t,dt,L,M)
    return M*(u_tp1-u_t)/dt+0.5*L*(u_tp1+u_t)
end

function P_J_system_Jm1(v_t,dt,L,M)
    coef_matrix=M/dt+0.5*L;
    rhs_J0=(M/dt-0.5*L)*v_t;
    return coef_matrix,rhs_J0
end

"""
    prolongation_graph(J::Int ,AD_EV::SparseMatrixCSC{Bool, Int64}, n::Int, m::Int)

Assemble prolongation matrix 'P' for prolongation from level J-1 to level 'J'.

"""
function prolongation_graph(J::Int, AD_EV::SparseMatrixCSC{Bool, Int64}, n::Int, m::Int)
    if J==1
        P_VV=I(n);
        #
        P_EV=0.5*AD_EV';
        # 
        P=vcat(P_VV,P_EV);
    else
        # P_EE
        P_EE=spzeros(0,0)
        P_e= sparse([if k in [2i-1, 2i+1] 0.5
                            elseif k ==2i 1.0 
                            else 0.0 end 
                            for i in 1:1:2.0^(J-1)-1, k in 1:1:2.0^J-1])'
        for _ = 1:m
            P_EE=blockdiag(P_EE,sparse(P_e))
        end
        # P_VV
        P_VV=I(n);
        # P_VE
        P_VE=zeros(n,size(P_EE)[2]);
        # P_EV
        P_EV=0.5*AD_EV';
        # 
        P=vcat(hcat(P_VV,P_VE),hcat(P_EV,P_EE));
    end
    return P
end

"""
    cn_mgm_iter!(u_iter::Matrix, J::Int, dt::Number, theta::Int, FE_matrix::SparseMatrixCSC, FE_mass_matrix::SparseMatrixCSC, n::Int, m::Int, rhs::Matrix, u_start::Matrix, nu1::Int, nu2::Int, mu::Int)

Perform Crank-Nicolson-Multigrid iteration to compute 'theta' time steps on 'u_iter' with step size 'dt'.

"""
function cn_mgm_iter!(u_iter::Matrix, J::Int, dt::Number, theta::Int, FE_matrix::SparseMatrixCSC, FE_mass_matrix::SparseMatrixCSC, n::Int, m::Int, rhs::Matrix, u_start::Matrix, nu1::Int, nu2::Int, mu::Int)

    ## smooth 
    coef_matrix_smooth = 0.5*FE_matrix+(1/dt)*FE_mass_matrix
    smooth_mats = SmoothingMatrices(Diagonal(1 ./diag(coef_matrix_smooth)), sparse(tril(coef_matrix_smooth,-1)), sparse(triu(coef_matrix_smooth,1)))

    for tau = 2:theta+1 # tau = 1 is current u_t, tau=k+1 is t+k*dt
        for _ = 1:nu1 
            u_iter_next = smooth(u_start[:,tau],(FE_mass_matrix/dt-0.5*FE_matrix)*u_iter[:,tau-1]+rhs[:,tau], smooth_mats)
            u_iter[:,tau] = copy(u_iter_next)
            u_start[:,tau] = copy(u_iter_next)
        end
    end

    d_tau = zeros(length(u_iter[:,1]),theta+1); # init
    for tau = 2:theta+1
        d_tau[:,tau] = P_J(u_iter[:,tau],u_iter[:,tau-1], dt, FE_matrix, FE_mass_matrix)-rhs[:,tau];
    end

    ## restriction operator

    AD_EV = ((FE_matrix[1:n,n+1:end] .!= 0));
    r = prolongation_graph(J, AD_EV, n, m)';

    ## restrict defects
    d_restr = zeros(length(r*d_tau[:,1]),theta+1)

    for tau = 1:theta+1
        d_restr[:,tau] = r*d_tau[:,tau]
    end

    ## restrict matrices
    FE_mass_matrix_Jm1 = r*FE_mass_matrix*r'
    FE_matrix_Jm1 = r*FE_matrix*r'

    ### coarse grid equations
    v_iter = zeros(length(r*d_tau[:,1]),theta+1)
        
    ## solve
    if J-1==0
        for tau=2:theta+1
            coef_matrix,rhs_J0 = P_J_system_Jm1(v_iter[:,tau-1],dt,FE_matrix_Jm1,FE_mass_matrix_Jm1)
            v_iter[:,tau] = coef_matrix\(d_restr[:,tau]+rhs_J0)
        end
    else
        for i=1:mu
            v_iter = cn_mgm_iter!(v_iter,J-1,dt,theta,FE_matrix_Jm1,FE_mass_matrix_Jm1,n,m,d_restr,v_iter,nu1,nu2,mu) ## ToDo: Fixe rechte Seite hinzufügen
        end
    end

    ### prolongation

    v_prolong = zeros(length(u_iter[:,1]),theta)
    for tau = 1:theta
        v_prolong[:,tau] = r'*v_iter[:,tau+1]
    end

    ### correction

    for tau = 2:theta+1
        u_iter[:,tau] = u_iter[:,tau]-v_prolong[:,tau-1]
    end

    u_iter_start = copy(u_iter)
    
    for tau = 2:theta+1 # tau = 1 is current u_t, tau=k+1 is t+k*dt
        for j = 1:nu2
            u_iter_next = smooth(u_iter[:,tau], (FE_mass_matrix/dt-0.5*FE_matrix)*u_iter[:,tau-1]+rhs[:,tau], smooth_mats)
            u_iter[:,tau] = copy(u_iter_next)
            u_iter_start[:,tau] = copy(u_iter_next)
        end
    end

    return u_iter
end

"""
    cn_mgm_iter!(u_iter::Matrix, J::Int, dt::Number, FE_matrix::SparseMatrixCSC, FE_mass_matrix::SparseMatrixCSC, n::Int, m::Int, rhs::Matrix, u_start::Matrix, nu1::Int, nu2::Int, mu::Int)

Perform Crank-Nicolson-Multigrid iteration to compute 'theta' time steps on 'u_iter' with step size 'dt'.

"""
function cn_mgm_iter!(u_iter::Vector, J::Int, dt::Number, FE_matrix::SparseMatrixCSC, FE_mass_matrix::SparseMatrixCSC, n::Int, m::Int, rhs::Vector, u_start::Vector, nu1::Int, nu2::Int, mu::Int)

    ## smooth 
    coef_matrix_smooth = 0.5*FE_matrix+(1/dt)*FE_mass_matrix
    smooth_mats = SmoothingMatrices(Diagonal(1 ./diag(coef_matrix_smooth)), sparse(tril(coef_matrix_smooth,-1)), sparse(triu(coef_matrix_smooth,1)))

    for _ = 1:nu1 
        u_iter_next = smooth(u_start,(FE_mass_matrix/dt-0.5*FE_matrix)*u_iter+rhs, smooth_mats)
        u_iter = copy(u_iter_next)
        u_start = copy(u_iter_next)
    end

    d_tau = zeros(length(u_iter)); # init
    d_tau = P_J(u_iter,u_iter, dt, FE_matrix, FE_mass_matrix)-rhs;

    ## restriction operator

    AD_EV = ((FE_matrix[1:n,n+1:end] .!= 0));
    r = prolongation_graph(J, AD_EV, n, m)';

    ## restrict defects
    d_restr = zeros(length(r*d_tau))

    d_restr = r*d_tau

    ## restrict matrices
    FE_mass_matrix_Jm1 = r*FE_mass_matrix*r'
    FE_matrix_Jm1 = r*FE_matrix*r'

    ### coarse grid equations
    v_iter = zeros(length(r*d_tau))
        
    ## solve
    if J-1==0
        coef_matrix,rhs_J0 = P_J_system_Jm1(v_iter,dt,FE_matrix_Jm1,FE_mass_matrix_Jm1)
        v_iter = coef_matrix\(d_restr+rhs_J0)
    else
        for i=1:mu
            v_iter = cn_mgm_iter!(v_iter,J-1,dt,FE_matrix_Jm1,FE_mass_matrix_Jm1,n,m,d_restr,v_iter,nu1,nu2,mu) ## ToDo: Fixe rechte Seite hinzufügen
        end
    end

    ### prolongation

    v_prolong = zeros(length(u_iter))
    v_prolong = r'*v_iter

    ### correction

    u_iter = u_iter-v_prolong

    u_iter_start = copy(u_iter)
    
    for j = 1:nu2
        u_iter_next = smooth(u_iter, (FE_mass_matrix/dt-0.5*FE_matrix)*u_iter+rhs, smooth_mats)
        u_iter = copy(u_iter_next)
        u_iter_start = copy(u_iter_next)
    end

    return u_iter
end

