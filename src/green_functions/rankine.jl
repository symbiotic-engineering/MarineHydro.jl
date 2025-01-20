struct Rankine <: GreensFunction end

_distance(x̅::PointNonDim, ξ̅::PointNonDim) = √∑((x̅ .- ξ̅) .^ 2)

# All functions below depends on wavenumber but the variable is actually unused.
# It is just part of the signature to have the same signature than the wave part of the Green function.

function greens(::Rankine, element_1, element_2, wavenumber=nothing)
    r = _distance(element_1[:center], element_2[:center])
    return 1 / r
end

function gradient_greens(::Rankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)

    r = _distance(element_1[:center], element_2[:center])
    r̅ = element_2[:center] - element_1[:center]
    ∇ₓ = 1/r^3*r̅
    if with_respect_to_first_variable
        return ∇ₓ
    else  # Gradient with respect to second variable
        return -∇ₓ
    end
end

function integral(::Rankine, element_1, element_2, wavenumber=nothing)
    @views begin
    point = element_1[:center]
    source_point = element_2[:center]
    source_vertices = element_2[:vertices]
    source_radius = element_2[:radius]
    source_normal = element_2[:normal]
    source_area = element_2[:area]
    end
    r = _distance(point, source_point)
    if r > 7 * source_radius # if far -> Rankine Direct
        function_integral = source_area / r
    else  # else (if close) deal with singularity
        function_integral = zero(source_area)
        r̄ = point - source_point
        GZ = dot(r̄, source_normal)
        RR = @inbounds [_distance(point, source_vertices[i,:]) for i in 1:4]
        NEXT_NODE = [2, 3, 4, 1]
        @views begin
        for L in 1:4
            DK = _distance(source_vertices[NEXT_NODE[L],:], source_vertices[L,:])
            if DK >= 1e-3 * source_radius
                PJ = (source_vertices[NEXT_NODE[L],:] - source_vertices[L,:]) / DK
                # The following GYX are called (a,b,c) in [Del]
                GYX = cross(source_normal, PJ)
                GY = dot(point - source_vertices[L,:], GYX)
                ANT = 2 * GY * DK
                DNT = (RR[NEXT_NODE[L]] + RR[L])^2 - DK^2 +
                    2 * abs(GZ) * (RR[NEXT_NODE[L]] + RR[L])
                ANL = RR[NEXT_NODE[L]] + RR[L] + DK
                DNL = RR[NEXT_NODE[L]] + RR[L] - DK
                ALDEN = log(ANL / DNL)
                if abs(GZ) >= 1e-4 * source_radius
                    AT = atan(ANT, DNT)  #check difference
                else
                    AT = zero(ANT) # zero of type 
                end
                if abs(GY) < 1e-5 # edge case
                    function_integral -= 2 * AT * abs(GZ)
                else
                    function_integral += GY * ALDEN - 2 * AT * abs(GZ)
                end
            end
        end
    end
end
    return function_integral
end


function integral_gradient(::Rankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    @views begin
    point = element_1[:center]
    source_point = element_2[:center]
    source_vertices = element_2[:vertices]
    source_normal = element_2[:normal]
    source_radius = element_2[:radius]
    source_area = element_2[:area]
        end
    r = _distance(point, source_point)
    r̅ =  point - source_point
    if r > 7 * source_radius # if far -> Rankine Direct
        derivative = 1 / r^3 * r̅
        derivative_integral = derivative * source_area
    else  # else (if close) deal with singularity
        r̄ = point - source_point
        GZ = dot(r̄, source_normal)
        RR = @inbounds [_distance(point, source_vertices[i, :]) for i in 1:4]
        DRX = @inbounds [(point - source_vertices[i, :]) / RR[i] for i in 1:4]
        
        VS_contribs = @inbounds [
            begin
                next_node = L % 4 + 1
                DK = _distance(source_vertices[next_node, :], source_vertices[L, :])
        
                if DK >= 1e-3 * source_radius
                    PJ = (source_vertices[next_node, :] - source_vertices[L, :]) / DK
        
                    # The following GYX are called (a, b, c) in [Del]
                    GYX = cross(source_normal, PJ)
                    GY = dot(point - source_vertices[L, :], GYX)
        
                    ANT = 2 * GY * DK
                    DNT = (RR[next_node] + RR[L])^2 - DK^2 +
                        2 * abs(GZ) * (RR[next_node] + RR[L])
                    ANL = RR[next_node] + RR[L] + DK
                    DNL = RR[next_node] + RR[L] - DK
                    ALDEN = log(ANL / DNL)
        
                    if abs(GZ) >= 1e-4 * source_radius
                        AT = atan(ANT, DNT)
                    else
                        AT = zero(ANT)
                    end
        
                    ANLX = DRX[next_node] + DRX[L]
                    ANTX = 2 * DK * GYX
                    DNTX = 2 * (RR[next_node] + RR[L] + abs(GZ)) * ANLX +
                        2 * sign(GZ) * (RR[next_node] + RR[L]) * source_normal
        
                    ALDEN .* GYX .- 2 * sign(GZ) * AT .* source_normal .+
                        GY * (DNL - ANL) ./ (ANL .* DNL) .* ANLX .-
                        2 * abs(GZ) .* (ANTX .* DNT .- DNTX .* ANT) ./ (ANT .* ANT .+ DNT .* DNT)
                else
                    zeros(typeof(source_area), 3)  # No contribution
                end
            end
            for L in 1:4
        ]
            # Aggregate contributions without mutation
            derivative_integral = -sum(VS_contribs)
            end

        if with_respect_to_first_variable
            return -derivative_integral
        else
            return derivative_integral
        end
    end
