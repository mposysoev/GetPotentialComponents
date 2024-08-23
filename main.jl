using ArgParse

include("src/GetPotentialComponents.jl")

using .GetPotentialComponents

# Define G2 functions list as a constant
const G2_FUNCTIONS_LIST = [G2(0.125, 7.0, 0.00),
    G2(4.000, 7.0, 3.00),
    G2(4.000, 7.0, 3.50),
    G2(4.000, 7.0, 4.00),
    G2(4.000, 7.0, 4.50),
    G2(4.000, 7.0, 5.00),
    G2(4.000, 7.0, 5.50),
    G2(4.000, 7.0, 6.00),
    G2(4.000, 7.0, 6.50)]

function main()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--component_type", "-c"
        help = "Type of component to calculate: 'pair', '3body' of 'diff'"
        arg_type = String
        required = true

        "--model_file", "-m"
        help = "Path to the model file"
        arg_type = String
        required = true

        "--output_file", "-o"
        help = "Path to the file for saving results"
        arg_type = String
        required = true

        "--r_max", "-r"
        help = "Maximum distance (radius) for calculations"
        arg_type = Float64
        required = true

        "--step_size", "-s"
        help = "Step size for distance increment in calculations"
        arg_type = Float64
        required = true

        "--separation", "-p"
        help = "Distance between two main atoms for 3-body component"
        arg_type = Float64
        required = false
    end

    parsed_args = parse_args(s)

    component_type = parsed_args["component_type"]
    model_file = parsed_args["model_file"]
    output_file = parsed_args["output_file"]
    r_max = parsed_args["r_max"]
    step_size = parsed_args["step_size"]
    separation = get(parsed_args, "separation", nothing)

    if component_type == "pair"
        println("Calculating pair component...")
        calculate_pair_component(
            model_file, output_file, r_max, step_size, G2_FUNCTIONS_LIST)
        println("Results saved in $output_file")
    elseif component_type == "3body"
        if separation === nothing
            error("The 'separation' parameter is required for calculating the 3-body component.")
        end
        println("Calculating 3-body component...")
        calculate_3body_component(
            model_file, output_file, r_max, step_size, separation, G2_FUNCTIONS_LIST)
        println("Results saved in $output_file")
    elseif component_type == "diff"
        if separation === nothing
            error("The 'separation' parameter is required for calculating the 3-body component.")
        end
        println("Calculating Difference of components...")
        calculate_3body_minus_2body_component(
            model_file, output_file, r_max, step_size, separation, G2_FUNCTIONS_LIST)
        println("Results saved in $output_file")
    else
        error("Invalid component type. Please use 'pair', '3body' of 'diff'.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
