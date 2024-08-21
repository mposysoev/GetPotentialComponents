using ArgParse

include("src/GetPotentialComponents.jl")

function main()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--component_type", "-c"
        help = "Type of component to calculate: 'pair' or '3body'"
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
        GetPotentialComponents.calculate_pair_component(
            model_file, output_file, r_max, step_size)
        println("Results saved to $output_file")
    elseif component_type == "3body"
        if separation === nothing
            error("Separation parameter is required for 3-body component calculation.")
        end
        println("Calculating 3-body component...")
        GetPotentialComponents.calculate_3body_component(
            model_file, output_file, r_max, step_size, separation)
        println("Results saved to $output_file")
    else
        error("Invalid component type. Please use 'pair' or '3body'.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
