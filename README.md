# Get Potential Components

This project aims to calculate pair and 3-body components of intermolecular potential from neural networks ([Flux.jl](https://fluxml.ai/)) that were trained on G2 symmetry functions.

## Installation

Firstly, you need to install all dependencies required for this project. You can do it easily just run:

```bash
make install
```

This command creates file `Manifest.toml`. It means that you are ready to run this code.

## Usage

There are two ways to run this code:

1. From `Jupyter Notebook` (recommended). For this, you need to run all cells from `run_me.ipynb` file. It will do everything for `example-methanol-model.bson` file automatically.
2. From the terminal. If you want to run this program as a command-line application:
   - `julia main.jl` – to see help
   - or `make run` to run the example for the `pair` component

## Authors and License

- Maksim Posysoev (<maxim.posysoev@gmail.com>)
- Prof. Alexander Lyubartsev (<alexander.lyubartsev@mmk.su.se>)

*GNU GPL-3.0 license*
