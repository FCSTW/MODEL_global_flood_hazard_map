import os
import platform
import json

def init():

    global CONFIG

    CONFIG = json.load(open('config.json', 'r'))

    # ================================================================================
    #
    # Modify the path
    #
    # ================================================================================

    # Get the root path of the project
    root_path                 = os.path.abspath(os.path.join(__file__, '../../../../'))
    CONFIG['path']['project'] = root_path

    # Loop over the path
    for key, value  in CONFIG['path'].items():

        # Replace {variable} in the string
        CONFIG['path'][key] = value.format(project=root_path)

        # Modify the path format
        if (platform.system() == 'Windows'):
            
            # If the path is starting with "/mnt/", convert it
            if (value[:5] == '/mnt/'):

                CONFIG['path'][key] = value.replace('/mnt/', '')

                # Convert the first slash to "://"
                CONFIG['path'][key] = CONFIG['path'][key].replace('/', '://', 1)

    return