# Docker and Singularity Container Setup

This folder contains configuration files and scripts for building Docker and Singularity containers for the project.
The containers encapsulate the application environment, in this case, the mpich runtime for EMOD and the 
fpg_observational_model is pre-installed with all necessary dependencies(Python 3.9, tskit, etc.). 

## Overview

**Docker** is a containerization platform that packages applications and their dependencies into isolated containers. 
It's ideal for local development and testing.

**Singularity** is a container platform designed for high-performance computing (HPC) environments. It's commonly used 
on Slurm clusters and our in-house Comps cluster, and provides better integration with HPC resources.

## Quick Start

**For Comps users:** Skip to [Using Comps](#1-using-comps-recommended-for-comps-users) section to build and push 
directly on Comps.

**For Slurm users:** Build a Singularity image following [Building Images](#building-images), then see 
[Running on Slurm](#2-running-on-slurm-clusters).

**For local development:** Build a Docker image following [Building Images](#building-images), then see 
[Local Development](#3-local-development-with-docker). The docker image can be used for ContainerPlatform as well. 

## Building Images

> [!IMPORTANT]
> If you are only running on Comps, you can skip this section and directly use the provided Singularity definition file
> to build the Singularity image on Comps. Please refer to the [Usage Scenarios](#usage-scenarios) section below for 
> more details.

### Docker Image

To build a Docker image from the Dockerfile:

> [!NOTE]
> You will need Docker installed. Please refer to the 
> [Docker Documentation](https://docs.docker.com/get-docker/) for installation instructions.

```bash
docker build -t <image-name>:<tag> -f <dockerfile-path> .
```

**Example:** Run from the root directory of the repository to build the Docker image using the provided Dockerfile:
```bash
docker build -t fpg-ob:latest -f docker/Dockerfile .
```

### Singularity Image

Singularity images can be built from Docker images or from Singularity definition files. After building, you will get 
a `.sif` file that can be used on HPC clusters. In this section, we focus on building Singularity images from 
definition files.

To build a Singularity image from a definition file (.def):

> [!NOTE]
> You will need Singularity (Apptainer) installed. Please refer to the 
> [Apptainer Documentation](https://apptainer.org/docs/user/main/quick_start.html) for installation 
> instructions. If you are running on Windows, consider using an Apptainer-in-Docker image to build Singularity images.
> Detailed Singularity build instructions on Windows can be found in the examples below.

```bash
apptainer build <image-name>.sif <definition-file>.def
```

**Examples:**

1. Using a Linux machine with Singularity/Apptainer installed, navigate to the repo root directory and run:
```bash
apptainer build ObsModel_rocky.sif docker/Singularity.def
```
Replace `apptainer` with `singularity` if you are using an older version of Singularity.

2. On Windows, you can use an Apptainer-in-Docker image to build the Singularity image, assuming you are in the
repo root directory:
```bash
docker run --rm `
  --device /dev/fuse `
  --security-opt seccomp=unconfined `
  --security-opt systempaths=unconfined `
  -v /var/run/docker.sock:/var/run/docker.sock `
  -v "${PWD}:/workspace" `
  ghcr.io/apptainer/apptainer:latest `
  apptainer build /workspace/ObsModel_rocky.sif /workspace/docker/Singularity.def
```
In this command:
- `--rm`: Automatically remove the container when it exits.
- `--device /dev/fuse`: Grants the container access to FUSE devices, necessary for Singularity operations.
- `--security-opt seccomp=unconfined` and `--security-opt systempaths=unconfined`: Relax security restrictions to allow 
Singularity to function properly within the Docker container.
- `-v /var/run/docker.sock:/var/run/docker.sock`: Mounts the Docker socket to allow the container to communicate with 
the Docker daemon on the host.
- `-v "${PWD}:/workspace"`: Mounts the current working directory to `/workspace` inside the container, allowing access 
to the definition file and output location.
- `ghcr.io/apptainer/apptainer:latest`: Specifies the Apptainer Docker image to use for building the Singularity image.
- `apptainer build /workspace/ObsModel_rocky.sif /workspace/docker/Singularity.def`: The command to build the Singularity image, 
specifying the output path and the definition file path.

In both cases, you will get a `ObsModel_rocky.sif` file in the current directory after the build is complete.

## Usage Scenarios

### 1. Using Comps (Recommended for Comps Users)

If you have access to Comps, you can push your Singularity image as an asset for easy distribution and version control.
> [!NOTE]
> You need to install idmtools and idmtools_platform_comps packages in your Python environment and have Comps access
> to use the script below. Run the following command to install the required packages if you haven't done so:
> ```bash
> pip install idmtools idmtools_platform_comps --index-url=https://packages.idmod.org/api/pypi/pypi-production/simple
> ```

Run the `push_singularity.py` script:

```bash
python push_singularity.py -f <path-to-image>
```

**Supported formats:**
- `.def` files (Singularity definition files) - builds the image on Comps
- `.sif` files (Singularity Image Format) - uploads pre-built image to Comps

**Optional parameters:**
```bash
python push_singularity.py -f <path-to-image> \
  --comps_url https://comps.idmod.org \
  --comps_env Calculon \
  --os_name rocky
```

> [!NOTE]
> Building from `.def` files on Comps is recommended for simplicity - you don't need to build locally or transfer 
> large `.sif` files. However, if you already have a tested `.sif` file, you can upload it directly for faster deployment.

> [!IMPORTANT]
> After successful upload, an asset ID file (e.g., `ObsModel_rocky.id`) will be created. Save this file as it 
> contains the asset ID needed to reference your image in Comps jobs. You can use this asset ID in your emodpy scripts 
> to run simulations and observational model in the Singularity image on Comps.
> A pre-built Singularity image is also available on Comps with asset ID stored in ['ObsModel_rocky.id'](ObsModel_rocky.id).

### 2. Running on Slurm Clusters

For Slurm-based HPC environments without Comps access:

1. Build your Singularity image locally or on a build node (see [Building Images](#building-images) section above)
2. Transfer the `.sif` file to your cluster if needed
3. Use the image in your Slurm job scripts.

### 3. Local Development with Docker

For local development, run a Docker container with the source code mounted. Assuming the image is built as `fpg-ob:latest` 
and you are in the `docker` folder:

```bash
docker run -it --rm `
  --mount type=bind,src="${PWD}/..",dst=/app `
  -w /app `
  fpg-ob:latest `
  bash
```

This command:
- Mounts the parent directory (repository root) to `/app` in the container
- Sets `/app` as the working directory
- Starts an interactive bash session

The `fpg_observational_model` package is pre-installed in the container. If you wish to develop and test changes to the 
package, you can reinstall it in editable mode.

Inside the container, install your project in development mode:

```bash
pip uninstall fpg_observational_model -y
pip install -e .
```

This allows you to edit code on your host machine while running it inside the container, with changes reflected 
immediately.

Another option is to use the `fpg-ob:latest` image to run with 
[ContainerPlatform](https://docs.idmod.org/projects/idmtools/en/latest/platforms/container/index.html)
for local testing of EMOD simulations with the observational model. Please refer to the ContainerPlatform documentation 
for details on how to set this up.

## Best Practices

- **Version your images**: Use meaningful tags (e.g., `v1.0`, `latest`, `dev`)
- **Keep containers lightweight**: Only include necessary dependencies
- **Document dependencies**: Maintain clear Dockerfile/definition files
- **Test locally first**: Validate your Docker image before building Singularity versions
- **Save asset IDs**: Keep the `.id` files generated by Comps uploads for future reference

## Troubleshooting

**Docker build fails**: Check that all required files are in the build context and paths are correct. Ensure you're 
running from the repository root directory.

**Singularity build requires sudo**: You may need root privileges or access to a build node with fakeroot support. On 
shared systems, contact your administrator.

**Permission issues with mounted volumes**: Ensure user IDs match between host and container. On Linux, you can use 
`--user $(id -u):$(id -g)` flag with Docker.

**Comps upload fails**: Verify your Comps credentials are configured correctly and you have network access to the Comps 
endpoint. Check that your `.def` or `.sif` file is valid and not corrupted.

**Asset ID file not created**: If `push_singularity.py` completes but no `.id` file appears, check the script output 
for errors during the Comps build process. The build may have failed on the Comps side.

**Windows Docker commands fail**: Ensure Docker Desktop is running and you're using PowerShell or Git Bash. Adjust path 
syntax if needed (`${PWD}` for PowerShell, `$(pwd)` for Git Bash).

## Additional Resources

- [Docker Documentation](https://docs.docker.com/)
- [Singularity Documentation](https://docs.sylabs.io/guides/latest/user-guide/)
- [Apptainer Documentation](https://apptainer.org/docs/)
