# © Anna Weller, University of Cologne, 2023

"""
    metric_tree_graph(; ℓ=1)

Create a tree graph with n=16 vertices and m=15 edges of lengths 'ℓ'.
"""
function metric_tree_graph(; ℓ=1)
    G_tree = path_graph(2);
    function add_children(G, v, r)
        n = nv(G)+1
        for _ in 1:r
            add_vertex!(G)
            add_edge!(G, v, n)
            n = n+1
        end
        return G
    end
    # first level
    G_tree = add_children(G_tree, 2, 2);
    # second level
    G_tree = add_children(G_tree, 3, 2);
    G_tree = add_children(G_tree, 4, 2);
    # third level
    G_tree = add_children(G_tree, 5, 2);
    G_tree = add_children(G_tree, 6, 2);
    G_tree = add_children(G_tree, 7, 2);
    G_tree = add_children(G_tree, 8, 2);
    #
    c_0 = 0.6;
    a = c_0/2;
    h_1 = sqrt(ℓ-c_0^2/4);
    c_1 = c_0+a;
    h_2 = sqrt(ℓ-c_1^2/4);
    c_2 = 2*a+2*c_0;
    h_3 = sqrt(ℓ-c_2^2/4);
    #
    coords=[
        [0, h_1+h_2+h_3+ℓ], #0
        [0, h_1+h_2+h_3],   #1
        [c_2/2, h_1+h_2],  #2 --> 4 & 5
        [-c_2/2, h_1+h_2], #3 --> 6 & 7
        [a/2+c_0/2, h_1],  #4 --> 8 & 9
        [a/2+c_0+a+c_0/2, h_1], #5 --> 10 & 11
        [-(a/2+c_0/2), h_1], #6
        [-(a/2+c_0+a+c_0/2), h_1], #7
        [a/2, 0], #8
        [a/2+c_0, 0], #9
        [a/2+c_0+a, 0], #10
        [a/2+c_0+a+c_0, 0], #11
        [-a/2, 0], #12
        [-(a/2+c_0), 0], #13
        [-(a/2+c_0+a), 0], #14
        [-(a/2+c_0+a+c_0), 0], #15
    ]
    return EquilateralMetricGraph(G_tree, ℓ, coords)
end

"""
    metric_cycle_graph(; ℓ=1)

Create a cycle graph with n=6 vertices and m=6 edges of lengths 'ℓ'.
"""
function metric_cycle_graph(; ℓ=1)
    G = CycleGraph(6);
    hoehe = ℓ/2*sqrt(3);
    coords = [
        [0, 2*ℓ],
        [h, 2*ℓ-hoehe],
        [hoehe, 2*ℓ-2*hoehe],
        [0, 0],
        [-hoehe, hoehe],
        [-hoehe, 2*hoehe],
    ];
    return EquilateralMetricGraph(G, ℓ, coords)
end

"""
    metric_graphene_graph(; ℓ=1)

Create a graphene graph with n=12 vertices and m=13 edges of lengths 'ℓ'.
"""
function metric_graphene_graph(; ℓ=1)
    G = CycleGraph(6);
    add_vertex!(G);
    add_edge!(G, 1, 7)
    for i in 8:12
        add_vertex!(G);
        add_edge!(G, i-1, i)
    end
    add_edge!(G, 12, 7)

    h = ℓ/2*sqrt(3);
    coords = [
        [0, 2*ℓ],
        [h, 2*ℓ-h],
        [h, 2*ℓ-2*h],
        [0, 0],
        [-h, h],
        [-h, 2*h],
        [0, 3*ℓ],
        [-h, 3*ℓ+h],
        [-h, 3*ℓ+2*h],
        [0, 3*ℓ+3*h],
        [h, 3*ℓ+2*h],
        [h, 3*ℓ+h],
    ];

    return EquilateralMetricGraph(G, ℓ, coords)
end

"""
    metric_star_graph(; ℓ=1)

Create a star graph with n=5 vertices and m=4 edges of lengths 'ℓ'.
"""
function metric_star_graph(; ℓ=1)
    coords = [[0, 0], [ℓ, 0], [-ℓ, 0], [0, ℓ], [0, -ℓ]]
    return EquilateralMetricGraph(star_graph(5), ℓ, coords)
end

"""
    metric_star_graph(ℓ_vec::Vector)

Create a star graph with edge lengths 'ℓ_vec'.
"""
function metric_star_graph(ℓ_vec::Vector)
    return MetricGraph(star_graph(length(ℓ_vec)+1), ℓ_vec, nothing)
end

"""
    metric_diamond_graph(; ℓ=1)

Create a diamond graph with n=4 vertices and m=5 edges of lengths 'ℓ'.
"""
function metric_diamond_graph(; ℓ=1)
    G = cycle_graph(4);
    add_edge!(G, 1, 3)
    hoehe = ℓ/2*sqrt(3)
    coords = [[-ℓ/2, 0], [0, hoehe], [ℓ/2, 0], [0, -hoehe]]
    Γ = EquilateralMetricGraph(G, ℓ, coords)
end

"""
    metric_lollipop_graph(n1::Int, n2::Int; ℓ=1)

Create a lollipop graph with clique of size 'n1' connected by an edge to a path of size 'n2', equilateral edge lengths 'ℓ'.
"""
function metric_lollipop_graph(; n1=3, n2=2, ℓ=1)
    G = lollipop_graph(n1, n2)
    if n1==3 && n2==2
        hoehe=ℓ/2*sqrt(3)
        coords=[[0, 0], [ℓ, 0], [ℓ/2, hoehe], [ℓ/2, hoehe+ℓ], [ℓ/2, hoehe+2*ℓ]]
        return EquilateralMetricGraph(G, ℓ, coords)
    elseif n1==3 && n2==1
        hoehe=ℓ/2*sqrt(3)
        coords=[[0, 0], [ℓ, 0], [ℓ/2, hoehe], [ℓ/2, hoehe+ℓ]]
        return EquilateralMetricGraph(G, ℓ, coords)
    else
        error(
            "metric_lollipop_graph is currently implemented for lollipop_graph(3,2) and lollipop_graph(3,1) only",
        )
    end
end

"""
    metric_barabasi_albert(n::Int, k::Int; ℓ=1, seed=nothing)

Create Barbási-Albert graph with 'n' vertices by growing an initial graph with 'k' vertices and
attaching each vertex with 'k' edges, see Graphs.barabasi_albert.

### Optional Arguments

  - ℓ: specify equilateral edge length, vector with edge length or ":non_equi" for random edge length
  - seed=nothing: set the RNG seed.
"""
function metric_barabasi_albert(n::Int, k::Int; ℓ=1, seed=nothing)
    G = barabasi_albert(n, k; seed=seed)
    if ℓ == :non_equi
        ℓ_vec = round.((rand(Uniform(1, 2), ne(G))); digits=2)
        return MetricGraph(G, ℓ_vec, nothing)
    elseif typeof(ℓ) == Vector
        return MetricGraph(G, ℓ_vec, nothing)
    else
        return EquilateralMetricGraph(G, ℓ, nothing)
    end
end

"""
    metric_erdos_renyi(n::Int, p::Number; ℓ=1)

Create equilateral Erdos-Renyi graph with 'n' vertices connected by edges with probability 'p', see Graphs.erdos_renyi.
"""
function metric_erdos_renyi(n::Int, p::Number; ℓ=1)
    G = erdos_renyi(n, p)
    return EquilateralMetricGraph(G, ℓ, nothing)
end

"""
    gaussian_initial(Γ; j_init = :randomly)

Assigne gaussian intial condition to one edge (randomly chosen or specified) of the graph.
"""
function gaussian_initial(Γ, j_init=:randomly)
    u0 = Vector{Function}(undef, ne(Γ.G))
    if j_init == :randomly
        j_init = rand(1:ne(Γ.G))
    end
    for j in 1:ne(Γ.G)
        if j == j_init
            x0 = Γ.ℓ/2;
            s = Γ.ℓ/10;
            u0[j] = x -> exp(-(x-x0)^2/s^2)
        else
            u0[j] = x -> 0
        end
    end
    return u0
end

"""
    model_initial_condition(Γ, j_init =:randomly)

Assigne intial condition with compact support to one edge (randomly chosen or specified) of the graph.
"""
function model_initial_condition(Γ, j_init=:randomly)
    u0 = Vector{Function}(undef, ne(Γ.G))
    if j_init == :randomly
        j_init = rand(1:ne(Γ.G))
    end
    for j in 1:ne(Γ.G)
        if j == j_init
            x0 = Γ.ℓ/2;
            s = Γ.ℓ/10;
            u0[j] = x -> exp(-(x-x0)^2/s^2)
        else
            u0[j] = x -> 0
        end
    end
    return u0
end
