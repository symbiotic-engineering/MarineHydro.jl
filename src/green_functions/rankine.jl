struct Rankine <: GreensFunction end

_distance(x, ξ) = norm(x .- ξ)

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
    point = element_1[:center]
    source_point = element_2[:center]
    source_vertices = element_2[:vertices]
    source_radius = element_2[:radius]
    source_normal = element_2[:normal]
    source_area = element_2[:area]
    r̄ = point - source_point
    r = norm(r̄)
    if r > 7 * source_radius # if far -> Rankine Direct
        integral = source_area / r
    else  # else (if close) deal with singularity
        integral = zero(source_area)
        GZ = dot(r̄, source_normal)
        RR = ntuple(i -> _distance(point, source_vertices[i,:]), 4)
        for index_vertex in 1:4
            current_vertex = source_vertices[index_vertex, :]
            index_next_vertex = mod1(index_vertex+1, 4)
            next_vertex = source_vertices[index_next_vertex, :]
            current_to_next = next_vertex - current_vertex
            segment_length = norm(current_to_next)
            if segment_length >= 1e-3 * source_radius
                GY = dot(point - current_vertex, cross(source_normal, current_to_next / segment_length))
                if abs(GZ) >= 1e-4 * source_radius
                    ANT = 2 * GY * segment_length
                    DNT = (RR[index_next_vertex] + RR[index_vertex])^2 - segment_length^2 +
                            2 * abs(GZ) * (RR[index_next_vertex] + RR[index_vertex])
                    integral -= 2*abs(GZ) * atan(ANT, DNT)
                end
                if abs(GY) > 1e-5
                    ANL = RR[index_next_vertex] + RR[index_vertex] + segment_length
                    DNL = RR[index_next_vertex] + RR[index_vertex] - segment_length
                    ALDEN = log(ANL / DNL)
                    integral += GY * ALDEN
                end
            end
        end
    end
    return integral
end

function integral_gradient(::Rankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    point = element_1[:center]
    source_point = element_2[:center]
    source_vertices = element_2[:vertices]
    source_normal = element_2[:normal]
    source_radius = element_2[:radius]
    source_area = element_2[:area]
    r̄ = point - source_point
    r = norm(r̄)
    if r > 7 * source_radius # if far -> Rankine Direct
        derivative = 1 / r^3 * r̄
        integral_gradient = derivative * source_area
    else  # else (if close) deal with singularity
        integral_gradient = zero(point)
        GZ = dot(r̄, source_normal)
        RR = ntuple(i -> _distance(point, source_vertices[i,:]), 4)
        DRX = ntuple(i -> (point - source_vertices[i,:]) / RR[i], 4)
        for index_vertex in 1:4
            current_vertex = source_vertices[index_vertex,:]
            index_next_vertex = mod1(index_vertex+1, 4)
            next_vertex = source_vertices[index_next_vertex, :]
            current_to_next = next_vertex - current_vertex
            segment_length = norm(current_to_next)
            if segment_length >= 1e-3 * source_radius
                GYX = cross(source_normal, current_to_next / segment_length)
                GY = dot(point - current_vertex, GYX)
                ANT = 2 * GY * segment_length
                DNT = (RR[index_next_vertex] + RR[index_vertex])^2 - segment_length^2 +
                    2 * abs(GZ) * (RR[index_next_vertex] + RR[index_vertex])
                ANL = RR[index_next_vertex] + RR[index_vertex] + segment_length
                DNL = RR[index_next_vertex] + RR[index_vertex] - segment_length
                ALDEN = log(ANL / DNL)
                if abs(GZ) >= 1e-4 * source_radius
                    AT = atan(ANT, DNT)  #check difference
                else
                    AT = zero(ANT)
                end
                #error bound error in the line below...
                ANLX = DRX[index_next_vertex] + DRX[index_vertex]
                ANTX = 2 * segment_length * GYX
                DNTX = 2 * (RR[index_next_vertex] + RR[index_vertex] + abs(GZ)) * ANLX +
                    2 * sign(GZ) * (RR[index_next_vertex] + RR[index_vertex]) * source_normal
                integral_gradient = integral_gradient .-
                                        (ALDEN .* GYX .- 2 * sign(GZ) * AT .* source_normal .+
                                        GY * (DNL - ANL) ./ (ANL .* DNL) .* ANLX .-
                                        2 * abs(GZ) .* (ANTX .* DNT .- DNTX .* ANT) ./ (ANT .* ANT .+ DNT .* DNT))
            end
        end
    end
    if with_respect_to_first_variable
        return -integral_gradient
    else  # Gradient with respect to second variable
        return integral_gradient
    end
end
