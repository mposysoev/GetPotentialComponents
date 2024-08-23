# Default target
.PHONY: all
all: run

# Install dependencies
.PHONY: install
install:
	julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run specific script (example)
.PHONY: run
run:
	julia --project=. main.jl -c pair -m example-methanol-model.bson -o pair-potential-component.txt -r 8 -s 0.1 -p 2.5

# Clean up any temporary or unwanted files
.PHONY: clean
clean:
	rm *.txt *.png

