module SeismicRayTrace
export raytrace
MAX_STEP = 10000
ϵ = 1.0e-5
α = 0.005

"""
set!(; maxit::Int=10000, epsilon::Float64=1.0e-5, alpha::Float64=0.1)

    set inner parameters
    maxit is maximum iteration step
    epsilon is stop condition of the iteration
    alpha is factor to reduce the step length of each iteration
"""
function set!(; maxit::Int=10000, epsilon::Float64=1.0e-5, alpha::Float64=0.005)
    global MAX_STEP = maxit
    global ϵ = epsilon
    global α = alpha
    return nothing
end

function _refraction_X(p::Float64, h::AbstractVector, v::AbstractVector)
    return p * sum(@.(v * h / sqrt(1.0 - p^2 * v^2)))
end

function _refraction_DpX(p::Float64, h::AbstractVector, v::AbstractVector)
    return sum(@.(v * h / (1.0 - p^2 * v^2)^(1.5)))
end

function _refraction_T(p::Float64, h::AbstractVector, v::AbstractVector)
    return sum(@.(h / (v * sqrt(1.0 - p^2 * v^2))))
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
    while abs(_refraction_X(p, h, v) - x0) > ϵ && step < MAX_STEP
        Δp = α * (x0 - _refraction_X(p, h, v)) / _refraction_DpX(p, h, v)
        p += min(1 / maxv - p - ϵ, max(-p + ϵ, Δp))
        step += 1
    end
    if step == MAX_STEP
        printstyled("MAX_STEP triggered, parameters:\nx0:", x0, "\nmodel_layer: ", join(string.(h), ' '),
                    "\nmodel_vel:", join(string.(v), ' '), "\n"; color=:yellow)
    end
    return p
end

"""
raytrace(depth1, depth2, Δ, model_depth, model_velocity, method)

    calculate travel time between two points in a layerd model

    depth1, depth2 is depth of two points
    Δ is horitonzal distance between these points
    model_depth and model_velocity should be vector
    method is a vector of `String`, in which can be any of refraction, reflection, guide
        if you want all type of phase, it can be reduced to ["all"]
"""
function raytrace(dep1::Float64, dep2::Float64, Δ::Float64, model_depth::Vector{Float64},
                  model_velocity::Vector{Float64}, method::Vector{T}) where {T<:AbstractString}
    # check parameter
    if findlast(<=(dep1), model_depth) == findlast(<=(dep2), model_depth)
        splitid = findlast(<=(dep1), model_depth)
        tdep = [model_depth[1] - 0.5; model_depth[1:splitid]; (dep1 + dep2) / 2.0; model_depth[splitid+1:end]]
        tvel = [0.0; model_velocity[1:splitid]; model_velocity[splitid]; model_velocity[splitid+1:end]]
    else
        splitid = length(model_velocity) + 3
        tdep = [model_depth[1] - 0.5; model_depth]
        tvel = [0.0; model_velocity]
    end
    d1 = min(dep1, dep2)
    d2 = max(dep1, dep2)
    i1 = findfirst(>(d1), tdep)
    i2 = findlast(<(d2), tdep)
    if dep1 == dep2
        i1 = i2 = splitid + 2
    end
    result = NamedTuple{(:x, :t, :p, :l, :type),Tuple{Float64,Float64,Float64,Int,String}}[]
    if "refraction" in method || "all" in method
        h = diff([d1; tdep[i1:i2]; d2])
        v = tvel[i1-1:i2]
        p0 = _refraction_raytrace(Δ, h, v)
        if isnan(p0)
            push!(result, (x=Δ, t=Δ / v[1], p=1 / v[1], l=0, type="refraction"))
        else
            push!(result,
                  (x=_refraction_X(p0, h, v), t=_refraction_T(p0, h, v), p=p0, l=0, type="refraction"))
        end
    end
    if "reflection" in method || "all" in method
        for il = 2:length(tdep)
            if tdep[il] == d1 || tdep[il] == d2 || il == splitid + 2
                continue
            end
            if il < i1
                h = [-diff([d1; tdep[i1-1:-1:il]]); diff([tdep[il:i2]; d2])]
                v = [tvel[i1-1:-1:il]; tvel[il:i2]]
            elseif il > i2
                h = [diff([d1; tdep[i1:il]]); -diff([tdep[il:-1:i2]; d2])]
                v = [tvel[i1-1:il-1]; tvel[il-1:-1:i2-1]]
            else
                continue
            end
            p0 = _refraction_raytrace(Δ, h, v)
            push!(result,
                  (x=_refraction_X(p0, h, v), t=_refraction_T(p0, h, v), p=p0,
                   l=il > splitid + 2 ? il - 2 : il - 1,
                   type="reflection"))
        end
    end
    if "guide" in method || "all" in method
        for il = 3:length(tdep)
            if il == splitid + 2
                continue
            end
            if il < i1
                if tvel[il-1] < tvel[il]
                    continue
                end
                p0 = 1.0 / tvel[il-1]
                if any(>=(tvel[il-1]), tvel[il:i2])
                    continue
                end
                h1 = diff([tdep[il:i1-1]; d1])
                v1 = tvel[il:i1-1]
                h2 = diff([tdep[il:i2-1]; d2])
                v2 = tvel[il:i2-1]
            elseif il > i2
                if tvel[il-1] > tvel[il]
                    continue
                end
                p0 = 1.0 / tvel[il]
                if any(>=(tvel[il]), tvel[i1-1:il-1])
                    continue
                end
                h1 = diff([d1; tdep[i1:il]])
                v1 = tvel[i1-1:il-1]
                h2 = diff([d2; tdep[i2:il]])
                v2 = tvel[i2-1:il-1]
            else
                continue
            end
            x1 = _refraction_X(p0, h1, v1)
            x2 = _refraction_X(p0, h2, v2)
            if Δ < (x1 + x2)
                continue
            end
            t1 = _refraction_T(p0, h1, v1)
            t2 = _refraction_T(p0, h2, v2)
            t3 = (Δ - x1 - x2) * p0
            push!(result, (x=Δ, t=t1 + t2 + t3, p=p0, l=il > splitid + 2 ? il - 2 : il - 1, type="guide"))
        end
    end
    return result
end

function raytrace(dep1::Real, dep2::Real, Δ::Real, model_depth::Vector{<:Real},
                  model_velocity::Vector{<:Real}, method::Union{AbstractString,Vector{<:AbstractString}})
    if typeof(method) <: AbstractString
        method = [method]
    end
    return raytrace(Float64(dep1), Float64(dep2), Float64(Δ), Float64.(model_depth),
                    Float64.(model_velocity), String.(method))
end
end
