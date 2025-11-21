from idmtools.core.platform_factory import Platform
from idmtools_platform_comps.utils.singularity_build \
                                         import SingularityBuildWorkItem
import os

def make_asset(def_file_name: str = 'Singularity.def', sif_file_name: str = None,
               comps_url: str = 'https://comps.idmod.org', comps_env: str = 'Calculon', os_name: str = 'rocky'
               ) -> None:
    """
    Build the EMOD observation model singularity image from definition file or singularity file on COMPS and save the asset id to a file.
    """
    # Prepare the platform
    plat_obj = Platform(type='COMPS',
                        endpoint=comps_url,
                        environment=comps_env)

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
    make_asset()
