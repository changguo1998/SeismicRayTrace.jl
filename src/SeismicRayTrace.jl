module SeismicRayTrace
export raytrace, raytrace_fastest, raytrace_guide
MAX_STEP = 100000
ϵ = 1.0e-5
α = 0.1

struct Phase
    layer::Vector{Int}
    polarization::Vector{Int}
end

"""
set!(; maxit::Int=100000, epsilon::Float64=1.0e-5, alpha::Float64=0.1)

    set inner parameters
    maxit is maximum iteration step
    epsilon is stop condition of the iteration
    alpha is factor to reduce the step length of each iteration
"""
function set!(; maxit::Int=NaN, epsilon::Float64=NaN, alpha::Float64=NaN)
    global MAX_STEP, ϵ, α
    if !isnan(maxit)
        MAX_STEP = maxit
    end
    if !isnan(epsilon)
        ϵ = epsilon
    end
    if !isnan(alpha)
        α = alpha
    end
    return nothing
end

function _refraction_X(p::Float64, h::AbstractVector, v::AbstractVector)
    return p * sum(@.(v * h / sqrt(1.0 - p^2 * v^2)))
end

function _refraction_DpX(p::Float64, h::AbstractVector, v::AbstractVector)
    return sum(@.(v * h / (1.0 - p^2 * v^2)^(1.5)))
end

function _refraction_D2pX(p::Float64, h::AbstractVector, v::AbstractVector)
    return 3.0 * p * sum(@.(v^3 * h / (1.0 - p^2 * v^2)^(2.5)))
end

function _refraction_T(p::Float64, h::AbstractVector, v::AbstractVector)
    return sum(@.(h / (v * sqrt(1.0 - p^2 * v^2))))
end

function _refraction_DpT(p::Float64, h::AbstractVector, v::AbstractVector)
    return p * sum(@.(v * h / (1.0 - p^2 * v^2)^(1.5)))
end

function _refraction_D2pT(p::Float64, h::AbstractVector, v::AbstractVector)
    return sum(@.(v * h * (1.0 + 2.0 * p^2 * v^2) / (1.0 - p^2 * v^2)^(2.5)))
end

# function _refraction_D2tX(p::Float64, h::AbstractVector, v::AbstractVector)
#     y′ = _refraction_DpT(p, h, v)
#     x′ = _refraction_DpX(p, h, v)
#     y′′ = _refraction_D2pT(p, h, v)
#     x′′ = _refraction_D2pX(p, h, v)
#     return (y′′ * x′ - x′′ * y′) / x′′^3
# end
function _refraction_Dp_TpoverXp(p::Float64, h::AbstractVector, v::AbstractVector)
    y′ = _refraction_DpT(p, h, v)
    x′ = _refraction_DpX(p, h, v)
    y′′ = _refraction_D2pT(p, h, v)
    x′′ = _refraction_D2pX(p, h, v)
    return (y′′ * x′ - x′′ * y′) / x′^2
end

function _refraction_raytrace(x0::Float64, h::AbstractVector, v::AbstractVector)
    if sum(h) == 0.0
        return NaN
    end
    if length(h) == 1
        return x0 / sqrt(x0^2 + h[1]^2) / v[1]
    end
    maxv = maximum(v)
    p = 0.0
    step = 1
    global α
    α0 = α
    while (abs(_refraction_X(p, h, v) - x0) > ϵ*x0) && (step < MAX_STEP)
        Δp = α0 * (x0 - _refraction_X(p, h, v)) / _refraction_DpX(p, h, v)
        if (Δp <= -p) || ((Δp + p) * maxv >= 1.0)
            α0 /= 2.0
        else
            p += Δp
            step += 1
            if α0 < α
                α0 *= 2.0
            end
        end
    end
    if step == MAX_STEP
        printstyled("MAX_STEP triggered, parameters:\nx0:", x0, "\nmodel_layer: ", join(string.(h), ' '),
                    "\nmvel:", join(string.(v), ' '), "\n"; color=:yellow)
    end
    return p
end

function _guide_raytrace(x0::Float64, v0::Float64, h::AbstractVector, v::AbstractVector)
    p0 = _refraction_raytrace(x0, h, v)
    if (_refraction_DpT(p0, h, v) / _refraction_DpX(p0, h, v)) <= (1.0 / v0)
        return NaN
    end
    p = p0
    step = 1
    global α
    α0 = α
    while abs(_refraction_DpT(p, h, v) / _refraction_DpX(p, h, v) - 1.0 / v0) > ϵ && step < MAX_STEP
        # println("p: ", p, ", X: ", _refraction_X(p, h, v), ", T: ", _refraction_T(p, h, v), ", DpT/DpX: ",
        #         _refraction_DpT(p, h, v) / _refraction_DpX(p, h, v))
        δp = α0 * (1.0 / v0 - _refraction_DpT(p, h, v) / _refraction_DpX(p, h, v)) / _refraction_Dp_TpoverXp(p, h, v)
        if (δp <= -p) || (δp >= (p0 - p))
            α0 /= 2.0
        else
            p += δp
            step += 1
            if α0 < α
                α0 *= 2.0
            end
        end
    end
    return p
end

function _guide_T(p::Float64, x0::Float64, v0::Float64, h::AbstractVector, v::AbstractVector)
    x = _refraction_X(p, h, v)
    t1 = _refraction_T(p, h, v)
    t2 = (x0 - x) / v0
    return t1 + t2
end

function _layerid(h::Vector{Float64}, h0::Float64, dep::Float64)
    c = h0
    i = 1
    while c <= dep
        i += 1
        c += h[i]
        if i == length(h)
            break
        end
    end
    return i
end

function _splitmodel(d1::Float64, d2::Float64, dep::Vector{Float64}, vel::VecOrMat{Float64})
    nvel = size(vel, 2)
    nl = length(dep)
    l1 = findlast(<=(d1), dep)
    l2 = findlast(<=(d2), dep)
    if l1 == l2
        g = vcat((dep[end] - dep[1]) * 0.1, diff([dep[1:l1]; d1; d2; dep[l1+1:end]]), 0.0)
        u = vcat(zeros(1, nvel), vel[[1:l1; l1; l1; l1+1:nl], :])
    else
        g = vcat((dep[end] - dep[1]) * 0.1, diff([dep[1:l1]; d1; dep[l1+1:l2]; d2; dep[l2+1:end]]), 0.0)
        u = vcat(zeros(1, nvel), vel[[1:l1; l1; l1+1:l2; l2; l2+1:nl], :])
    end

    nzerothick = count(iszero, g)
    h = zeros(length(g) - nzerothick + 1)
    v = zeros(length(g) - nzerothick + 1, nvel)
    p = 1
    for i = 1:length(g)-1
        if g[i] != 0
            h[p] = g[i]
            v[p, :] .= u[i, :]
            p += 1
        end
    end
    h[p] = g[end]
    v[p, :] .= u[end, :]

    i1 = _layerid(h, dep[1], d1)
    i2 = _layerid(h, dep[1], d2)
    newlayer = zeros(Int, length(dep))
    for i = 1:length(dep)
        newlayer[i] = _layerid(h, dep[1], dep[i])
    end
    return (h, v, i1, i2, newlayer)
end

function _eqmodel(h::Vector{Float64}, v::VecOrMat{Float64}, path::Vector{Int}, phase::Vector{Int})
    nlayer = 0
    for i = 1:length(path)-1
        nlayer += abs(path[i+1] - path[i])
    end
    l = zeros(nlayer)
    w = zeros(nlayer)
    p = 1
    for i = 1:length(path)-1
        l1 = path[i]
        l2 = path[i+1]
        ip = phase[i]
        if l1 < l2
            for j = l1:l2-1
                l[p] = h[j]
                w[p] = v[j, ip]
                p += 1
            end
        else
            for j = l1-1:-1:l2
                l[p] = h[j]
                w[p] = v[j, ip]
                p += 1
            end
        end
    end
    return (l, w)
end

function _new_path(path::Vector{Int}, table::Vector{Int}, i1::Int, i2::Int)
    newpath = zeros(Int, length(path) + 2)
    newpath[1] = i1
    newpath[end] = i2
    for i = 1:length(path)
        newpath[i+1] = table[path[i]]
    end
    return newpath
end

function _cal_time(h, v, x0)
    p = _refraction_raytrace(x0, h, v)
    t = _refraction_T(p, h, v)
    xt = _refraction_X(p, h, v)
    return (p=p, t=t, x=xt)
end

function _cal_guide(h, v, x0, v0)
    if any(>=(v0), v)
        return (p=NaN, t=NaN, x=NaN)
    end
    p = _guide_raytrace(x0, v0, h, v)
    if isnan(p)
        return (p=NaN, t=NaN, x=NaN)
    end
    t = _guide_T(p, x0, v0, h, v)
    xt = x0
    return (p=p, t=t, x=xt)
end

# * core
function raytrace(dep1::Float64, dep2::Float64, Δ::Float64, mdep::Vector{Float64}, mvel::VecOrMat{Float64},
                  path::Vector{Vector{Int}}, pol::Vector{Vector{Int}})
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    (h, v, i1, i2, newlayer) = _splitmodel(d1, d2, mdep, mvel)
    return map(path, pol) do pa, po
        np = _new_path(pa, newlayer, i1, i2)
        (l, w) = _eqmodel(h, v, np, po)
        _cal_time(l, w, Δ)
    end
end

function raytrace(dep1::Float64, dep2::Float64, Δ::Float64, mdep::Vector{Float64}, mvel::VecOrMat{Float64},
                  path::Vector{Int}, pol::Vector{Int})
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    (h, v, i1, i2, newlayer) = _splitmodel(d1, d2, mdep, mvel)
    np = _new_path(path, newlayer, i1, i2)
    (l, w) = _eqmodel(h, v, np, pol)
    return _cal_time(l, w, Δ)
end

function raytrace_guide(dep1::Float64, dep2::Float64, Δ::Float64, mdep::Vector{Float64}, mvel::Vector{Float64},
                        layer::Vector{Int})
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    (h, v, i1, i2, newlayer) = _splitmodel(d1, d2, mdep, mvel)
    return map(layer) do nl
        path = Int[nl]
        pol = Int[1, 1]
        np = _new_path(path, newlayer, i1, i2)
        (l, w) = _eqmodel(h, v, np, pol)
        v0 = (nl == 1) ? mvel[1] : max(mvel[nl-1], mvel[nl])
        _cal_guide(l, w, Δ, v0)
    end
end

function raytrace_guide(dep1::Float64, dep2::Float64, Δ::Float64, mdep::Vector{Float64}, mvel::Vector{Float64},
                        layer::Int)
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    (h, v, i1, i2, newlayer) = _splitmodel(d1, d2, mdep, mvel)
    path = Int[layer]
    pol = Int[1, 1]
    np = _new_path(path, newlayer, i1, i2)
    (l, w) = _eqmodel(h, v, np, pol)
    # println("l: ", l)
    # println("w: ", w)
    v0 = (layer == 1) ? mvel[1] : max(mvel[layer-1], mvel[layer])
    # println("v0: ", v0)
    return _cal_guide(l, w, Δ, v0)
end

# * adapter
function raytrace(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVecOrMat{<:Real},
                  path::AbstractVector{<:Integer}, pol::AbstractVector{<:Integer})
    return raytrace(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), Int.(path), Int.(pol))
end

function raytrace(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVecOrMat{<:Real},
                  path::AbstractVector{<:Integer})
    pol = ones(Int, length(path) + 1)
    return raytrace(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), Int.(path), pol)
end

function raytrace(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVecOrMat{<:Real},
                  path::AbstractVector{<:AbstractVector{<:Integer}}, pol::AbstractVector{<:AbstractVector{<:Integer}})
    npath = Vector{Int}[]
    npol = Vector{Int}[]
    for i in eachindex(path)
        push!(npath, Int.(path[i]))
        push!(npol, Int.(pol[i]))
    end
    return raytrace(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), npath, npol)
end

function raytrace(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVecOrMat{<:Real},
                  path::AbstractVector{<:AbstractVector{<:Integer}})
    npath = Vector{Int}[]
    npol = Vector{Int}[]
    for i in eachindex(path)
        push!(npath, Int.(path[i]))
        push!(npol, ones(Int, length(path[i]) + 1))
    end
    return raytrace(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), npath, npol)
end

function raytrace_guide(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVector{<:Real},
                        layer::AbstractVector{<:Integer})
    return raytrace_guide(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), Int.(layer))
end

function raytrace_guide(dep1::Real, dep2::Real, Δ::Real, mdep::AbstractVector{<:Real}, mvel::AbstractVector{<:Real},
                        layer::Integer)
    return raytrace_guide(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(mdep), Float64.(mvel), Int(layer))
end

function raytrace_fastest(dep1::Real, dep2::Real, Δ::Real, mdep::Vector{<:Real}, mvel::Vector{<:Real})
    tmin = Inf
    reflectlayer = Vector{Int}[]
    gdep = Int[]
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    push!(reflectlayer, Int[])
    for i in eachindex(mdep)
        if (mdep[i] < d1) || (mdep[i] > d2)
            push!(reflectlayer, Int[i])
            push!(gdep, i)
        end
    end
    phases = [raytrace(dep1, dep2, Δ, mdep, mvel, reflectlayer);
              raytrace_guide(dep1, dep2, Δ, mdep, mvel, gdep)]
    tvec = Float64[]
    idvec = Int[]
    for i in eachindex(phases)
        if !isnan(phases[i].t)
            push!(tvec, phases[i].t)
            push!(idvec, i)
        end
    end
    (_, mint_id) = findmin(tvec)
    id = idvec[mint_id]
    if id <= length(reflectlayer)
        ptype = "reflect"
        path = reflectlayer[id]
    else
        ptype = "guide"
        path = gdep[id-length(reflectlayer)]
    end
    return (phase = phases[id], path=path, type=ptype)
end

end
