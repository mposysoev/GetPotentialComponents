module GetPotentialComponents

using BSON: @load
using Flux

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

function main()
    calculate_pair_component("methanol-CG-NN.bson", "output_results.txt", 10.0, 0.1)
end

end # module GetPotentialComponents
