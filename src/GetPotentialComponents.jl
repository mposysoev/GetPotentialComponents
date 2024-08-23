module GetPotentialComponents

using BSON: @load
using Flux

export calculate_pair_component, calculate_3body_component,
       calculate_3body_minus_2body_component, main, G2

struct G2
    eta::Float64
    rcutoff::Float64
    rshift::Float64
end

function load_model(file_path::AbstractString)::Flux.Chain
    model = nothing
    @load file_path model
    return model
end

function init_positions_matrix(size::Int64)
    return zeros(size, 3)
end

function distance(positions, i, j)
    return (positions[i, 1] - positions[j, 1])^2 + (positions[i, 2] - positions[j, 2])^2 +
           (positions[i, 3] - positions[j, 3])^2
end

function init_distance_matrix(positions)
    N = size(positions, 1)
    dist_mat = zeros(Float64, (N, N))

    for i in 1:N
        for j in 1:N
            dist_mat[i, j] = distance(positions, i, j)
        end
    end

    return dist_mat
end

function distance_cutoff(distance::Float64, rcutoff = 6.0)::Float64
    if distance > rcutoff
        return 0.0
    else
        return (0.5 * (cos(π * distance / rcutoff) + 1.0))
    end
end

function calculate_g2_element(distance, eta, rcutoff, rshift)::Float64
    if distance > 0.0
        return exp(-eta * (distance - rshift)^2) * distance_cutoff(distance, rcutoff)
    else
        return 0.0
    end
end

function calculate_g2_function(distances, eta, rcutoff, rshift)::Float64
    sum = 0.0
    for distance in distances
        sum += calculate_g2_element(distance, eta, rcutoff, rshift)
    end
    return sum
end

function calculate_g2_matrix(distance_matrix, g2_func_params_list::Vector{G2})
    N = size(distance_matrix)[1]
    g2_matrix = zeros(Float64, (N, length(g2_func_params_list)))
    for i in 1:N
        distance_vector = distance_matrix[i, :]
        for (j, g2_function) in enumerate(g2_func_params_list)
            eta = g2_function.eta
            rcutoff = g2_function.rcutoff
            rshift = g2_function.rshift
            g2_matrix[i, j] = calculate_g2_function(distance_vector, eta, rcutoff, rshift)
        end
    end

    return g2_matrix
end

function atomic_energy(inputlayer, model)
    E::Float64 = model(inputlayer)[1]
    return (E)
end

function calculate_system_energy_vector(symmFuncMatrix, model)
    N = size(symmFuncMatrix)[1]
    E = zeros(Float64, N)
    for i in 1:N
        E[i] = atomic_energy(symmFuncMatrix[i, :], model)
    end
    return E
end

function calculate_pair_component(
        model_file::String, output_file::String, r_max::Float64, step_size::Float64,
        g2_func_params_list::Vector{G2})
    model = load_model(model_file)
    positions = init_positions_matrix(2)
    distance_matrix = init_distance_matrix(positions)

    num_iter = r_max / step_size

    open(output_file, "w") do io
        println(io, "# Pair potential for $model_file")
        println(io, "# r [Å],  U(r)")

        for i in 1:num_iter
            distance_matrix[1, 2] += step_size
            distance_matrix[2, 1] += step_size

            g2_matrix = calculate_g2_matrix(distance_matrix, g2_func_params_list)
            energy_vector = calculate_system_energy_vector(g2_matrix, model)

            total_energy = sum(energy_vector)
            dist = round(distance_matrix[1, 2], digits = 5)

            println(io, "$dist    $total_energy")
        end
    end
end

function calculate_3body_component(
        model_file::String, output_file::String, r_max::Float64, step_size::Float64, separation::Float64,
        g2_func_params_list::Vector{G2})
    model = load_model(model_file)
    positions = init_positions_matrix(3)

    num_iter = round(Int, 2 * r_max / step_size)
    initial_x = positions[3, 1] - r_max  # save the initial value of X

    open(output_file, "w") do io
        println(io,
            "# 3-body potential component for $model_file with distance between particles = $(separation) Å")
        println(io, "# Grid (X,Y) from (-$(r_max), +$(r_max)) to (+$(r_max), -$(r_max))")
        println(io, "# X, Y, U(r)")

        # scanning particle
        positions[3, 1] = initial_x  # set X coordinate
        positions[3, 2] += r_max  # shift Y coordinate

        # separate two main particles
        positions[1, 1] -= separation / 2.0
        positions[2, 1] += separation / 2.0

        for y in 1:num_iter
            for x in 1:num_iter
                distance_matrix = init_distance_matrix(positions)
                g2_matrix = calculate_g2_matrix(distance_matrix, g2_func_params_list)
                energy_vector = calculate_system_energy_vector(g2_matrix, model)
                total_energy = sum(energy_vector)
                println(io, "$(positions[3, 1]), $(positions[3, 2]), $total_energy")

                # change scanning particle position
                positions[3, 1] += step_size
            end
            # return X coordinate
            positions[3, 1] = initial_x  # restore the initial value of X

            # step for Y coordinate
            positions[3, 2] -= step_size
        end
    end
end

function calculate_pair_component_once(
        model::Flux.Chain, distance::Float64, g2_func_params_list::Vector{G2})
    positions = init_positions_matrix(2)

    positions[1, 1] += distance

    distance_matrix = init_distance_matrix(positions)

    g2_matrix = calculate_g2_matrix(distance_matrix, g2_func_params_list)
    energy_vector = calculate_system_energy_vector(g2_matrix, model)

    total_energy = sum(energy_vector)

    return total_energy
end

function calculate_3body_minus_2body_component(
        model_file::String, output_file::String, r_max::Float64, step_size::Float64, separation::Float64,
        g2_func_params_list::Vector{G2})
    model = load_model(model_file)
    postions_3_particles = init_positions_matrix(3)
    postions_2_particles_left = init_positions_matrix(2)
    postions_2_particles_right = init_positions_matrix(2)

    num_iter = round(Int, 2 * r_max / step_size)
    initial_x = postions_3_particles[3, 1] - r_max  # save the initial value of X

    open(output_file, "w") do io
        println(io,
            "# 3-body potential component for $model_file with distance between particles = $(separation) Å")
        println(io, "# Grid (X,Y) from (-$(r_max), +$(r_max)) to (+$(r_max), -$(r_max))")
        println(io, "# X, Y, U(r)")

        # scanning particle
        postions_3_particles[3, 1] = initial_x  # set X coordinate
        postions_3_particles[3, 2] += r_max  # shift Y coordinate

        postions_2_particles_left[2, 1] = initial_x  # set X coordinate
        postions_2_particles_left[2, 2] += r_max  # shift Y coordinate

        postions_2_particles_right[2, 1] = initial_x  # set X coordinate
        postions_2_particles_right[2, 2] += r_max  # shift Y coordinate

        # separate two main particles
        postions_3_particles[1, 1] -= separation / 2.0
        postions_3_particles[2, 1] += separation / 2.0

        postions_2_particles_left[1, 1] -= separation / 2.0
        postions_2_particles_right[1, 1] += separation / 2.0

        for y in 1:num_iter
            for x in 1:num_iter
                distance_matrix_3_particles = init_distance_matrix(postions_3_particles)
                distance_matrix_2_left = init_distance_matrix(postions_2_particles_left)
                distance_matrix_2_right = init_distance_matrix(postions_2_particles_right)

                g2_matrix_3_particles = calculate_g2_matrix(
                    distance_matrix_3_particles, g2_func_params_list)
                g2_matrix_2_left = calculate_g2_matrix(
                    distance_matrix_2_left, g2_func_params_list)
                g2_matrix_2_right = calculate_g2_matrix(
                    distance_matrix_2_right, g2_func_params_list)

                energy_vector_3_particles = calculate_system_energy_vector(
                    g2_matrix_3_particles, model)
                energy_vector_2_left = calculate_system_energy_vector(
                    g2_matrix_2_left, model)
                energy_vector_2_right = calculate_system_energy_vector(
                    g2_matrix_2_right, model)

                # Из полной энергии с трёмя частицами,
                # вычесть энергию системы из двух частиц (сканирующей и статичной левой),
                # вычесть энергию системы из двух частиц (сканирующей и статичной правой),
                # вычесть энергию взаимодействия между двумя частицами.
                total_energy = sum(energy_vector_3_particles) - sum(energy_vector_2_left) -
                               sum(energy_vector_2_right) - calculate_pair_component_once(
                    model, separation, g2_func_params_list)

                println(io,
                    "$(postions_3_particles[3, 1]), $(postions_3_particles[3, 2]), $total_energy")

                # change scanning particle position
                postions_3_particles[3, 1] += step_size
                postions_2_particles_left[2, 1] += step_size
                postions_2_particles_right[2, 1] += step_size
            end
            # return X coordinate
            postions_3_particles[3, 1] = initial_x
            postions_2_particles_left[2, 1] = initial_x
            postions_2_particles_right[2, 1] = initial_x

            # step for Y coordinate
            postions_3_particles[3, 2] -= step_size
            postions_2_particles_left[2, 2] -= step_size
            postions_2_particles_right[2, 2] -= step_size
        end
    end
end

function main(
        component_type::String, model_file::String, output_file::String, r_max::Float64,
        step_size::Float64, g2_func_params_list::Vector{G2}, separation::Union{
            Float64, Nothing} = nothing)
    if component_type == "pair"
        println("Calculating pair component...")
        calculate_pair_component(
            model_file, output_file, r_max, step_size, g2_func_params_list)
        println("Results saved in $output_file")
    elseif component_type == "3body"
        if separation === nothing
            error("The 'separation' parameter is required for calculating the 3-body component.")
        end
        println("Calculating 3-body component...")
        calculate_3body_component(
            model_file, output_file, r_max, step_size, separation, g2_func_params_list)
        println("Results saved in $output_file")
    elseif component_type == "diff"
        if separation === nothing
            error("The 'separation' parameter is required for calculating the 3-body component.")
        end
        println("Calculating Difference of components...")
        calculate_3body_minus_2body_component(
            model_file, output_file, r_max, step_size, separation, g2_func_params_list)
        println("Results saved in $output_file")
    else
        error("Invalid component type. Please use 'pair', '3body' of 'diff'.")
    end
end

end # module GetPotentialComponents
