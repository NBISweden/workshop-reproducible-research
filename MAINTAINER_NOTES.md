# Notes

## Tutorial environment package management

Multi-platform package environment management for each tutorial is handled using Pixi.

### For maintainers

To upgrade packages, run:

```bash
pixi upgrade
```

Then update the conda specification files:

```bash
pixi run update-conda-specfiles
```

Alternatively, a GitHub Action automatically updates spec files if the Pixi files are updated.

Once Pixi is more widely adopted, we may consider switching to using Pixi directly for environment creation.

### For students

Students create a tutorial environment by running:

```bash
conda create --name <env_name> --file env/conda-specfiles/<env_name>-<platform>_conda_spec.txt
```

where `<env_name>` is the desired name for the conda environment and `<platform>` is one of the following options:

- `linux-64`
- `osx-64` (MacOS Intel CPUs)
- `osx-arm64` (MacOS M-Series CPUs)
