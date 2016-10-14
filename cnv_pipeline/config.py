import os

import configparser


this_dir = os.path.dirname(os.path.realpath(__file__))
config = configparser.ConfigParser()


def load_config():
    """Read configuration file."""
    config.read(os.path.join(this_dir, 'config.ini'))
    return config