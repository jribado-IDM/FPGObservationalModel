import warnings

from idmtools.core.platform_factory import Platform
from idmtools_platform_comps.utils.singularity_build \
                                         import SingularityBuildWorkItem
import os

def make_asset(def_file_name: str = 'Singularity.def', sif_file_name: str = None,
               comps_url: str = 'https://comps.idmod.org', comps_env: str = 'Calculon', os_name: str = 'rocky'
               ) -> None:
    """
    Build the EMOD observation model singularity image from definition file or singularity file on COMPS and save the
    asset id to a file.
    Args:
        def_file_name (str): Path to the Singularity definition file(.def), default is 'Singularity.def'. If provided,
                             this will be used to build the image instead of sif_file_name.
        sif_file_name (str): Path to the Singularity image file(.sif), default is None which means build from def file.
        comps_url (str): COMPS endpoint URL, default is 'https://comps.idmod.org'.
        comps_env (str): COMPS environment name to use, default is 'Calculon'.
        os_name (str): Operating system name for the image, default is 'rocky'.
    """
    # Prepare the platform
    plat_obj = Platform(type='COMPS',
                        endpoint=comps_url,
                        environment=comps_env)

    if def_file_name is not None and sif_file_name is not None:
        warnings.warn("Both def_file_name and sif_file_name are provided. def_file_name will be used.")
        sif_file_name = None

    # Create Singularity build work item
    if def_file_name is not None:
        # Check if the definition file exists
        if not os.path.exists(def_file_name):
            raise FileNotFoundError(f'Specified SIF definition file {def_file_name} does not exist.')
        # Build from definition file
        sbwi_obj = SingularityBuildWorkItem(name='ObsModel_'+os_name,
                                            definition_file=def_file_name,
                                            image_name='ObsModel_'+os_name+'.sif',
                                            force=True)
    elif sif_file_name is not None:
        # Check if the singularity file exists
        if not os.path.exists(sif_file_name):
            raise FileNotFoundError(f'Specified SIF file {sif_file_name} does not exist.')
        # Build from existing singularity file
        sbwi_obj = SingularityBuildWorkItem(name='ObsModel_'+os_name,
                                            definition_content=sif_file_name,
                                            force=True)
    else:
        raise ValueError('Either def_file_name or sif_file_name must be provided.')

    # Wait until the image is built
    ac_obj = sbwi_obj.run(wait_until_done=True, platform=plat_obj)

    # Save asset id for sif to file
    ac_obj.to_id_file('ObsModel_'+os_name+'.id')

    return None


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build Observational Model Singularity image and save asset id to file.')
    parser.add_argument('--file', '-f', type=str, default='Singularity.def',
                        help='Path to the Singularity definition file (.def) or image file (.sif) to build from. '
                             'Default is Singularity.def')
    parser.add_argument('--comps_url', type=str, default='https://comps.idmod.org',
                        help='COMPS endpoint URL. Default is https://comps.idmod.org')
    parser.add_argument('--comps_env', type=str, default='Calculon',
                        help='COMPS environment name to use. Default is Calculon')
    parser.add_argument('--os_name', type=str, default='rocky',
                        help='Operating system name for the image. Default is rocky')
    args = parser.parse_args()
    if args.file.endswith('.def'):
        def_file = args.file
        sif_file = None
    elif args.file.endswith('.sif'):
        def_file = None
        sif_file = args.file
    else:
        raise ValueError('The provided file must be either a .def or .sif file.')
    make_asset(def_file_name=def_file, sif_file_name=sif_file, comps_url=args.comps_url,
               comps_env=args.comps_env, os_name=args.os_name)
