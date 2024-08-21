module GetPotentialComponents

using BSON: @load
using Flux

export calculate_pair_component, calculate_3body_component

struct G2
    eta::Float64
    rcutoff::Float64
    rshift::Float64
end

const G2_FUNCTIONS_LIST = [G2(0.125, 7.0, 0.00),
    G2(4.000, 7.0, 3.00),
    G2(4.000, 7.0, 3.50),
    G2(4.000, 7.0, 4.00),
    G2(4.000, 7.0, 4.50),
    G2(4.000, 7.0, 5.00),
    G2(4.000, 7.0, 5.50),
    G2(4.000, 7.0, 6.00),
    G2(4.000, 7.0, 6.50)]

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
        model_file::String, output_file::String, r_max::Float64, step_size::Float64)
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

            g2_matrix = calculate_g2_matrix(distance_matrix, G2_FUNCTIONS_LIST)
            energy_vector = calculate_system_energy_vector(g2_matrix, model)

            total_energy = sum(energy_vector)
            dist = round(distance_matrix[1, 2], digits = 5)

            println(io, "$dist    $total_energy")
        end
    end
end

function calculate_3body_component(
        model_file::String, output_file::String, r_max::Float64, step_size::Float64, separation::Float64)
    model = load_model(model_file)
    positions = init_positions_matrix(3)

    num_iter = round(Int, 2 * r_max / step_size)
    initial_x = positions[3, 1] - r_max  # сохранить начальное значение X

    open(output_file, "w") do io
        println(io,
            "# 3-body potential component for $model_file with distance between particles = $(separation) Å")
        println(io, "# Grid (X,Y) from (-$(r_max), +$(r_max)) to (+$(r_max), -$(r_max))")
        println(io, "# X, Y, U(r)")

        # scaning particle
        positions[3, 1] = initial_x  # установить X coordinate
        positions[3, 2] += r_max  # shift Y coordinate

        # separate two main particles
        positions[1, 1] -= separation / 2.0
        positions[2, 1] += separation / 2.0

        for y in 1:num_iter
            for x in 1:num_iter
                distance_matrix = init_distance_matrix(positions)
                g2_matrix = calculate_g2_matrix(distance_matrix, G2_FUNCTIONS_LIST)
                energy_vector = calculate_system_energy_vector(g2_matrix, model)
                total_energy = sum(energy_vector)
                println(io, "$(positions[3, 1]), $(positions[3, 2]), $total_energy")

                # change scaning particle position
                positions[3, 1] += step_size
            end
            # return X coordinate
            positions[3, 1] = initial_x  # восстановить начальное значение X

            # step for Y coordinate
            positions[3, 2] -= step_size
        end
    end
end

function main(
        component_type::String, model_file::String, output_file::String, r_max::Float64,
        step_size::Float64, separation::Union{Float64, Nothing} = nothing)
    if component_type == "pair"
        println("Вычисление парной компоненты...")
        calculate_pair_component(model_file, output_file, r_max, step_size)
        println("Результаты сохранены в $output_file")
    elseif component_type == "3body"
        if separation === nothing
            error("Для вычисления трехтелесной компоненты требуется параметр separation.")
        end
        println("Вычисление трехтелесной компоненты...")
        calculate_3body_component(model_file, output_file, r_max, step_size, separation)
        println("Результаты сохранены в $output_file")
    else
        error("Неверный тип компоненты. Пожалуйста, используйте 'pair' или '3body'.")
    end
end

end # module GetPotentialComponents
