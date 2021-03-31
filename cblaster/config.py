import configparser
import appdirs
import logging

from pathlib import Path


LOG = logging.getLogger(__name__)


def get_config_dir():
    path = appdirs.user_config_dir(appname="cblaster")
    return Path(path)


def get_config_parser():
    # Get configuration directory, platform agnostic
    config_dir = get_config_dir()
    if not config_dir.is_dir():
        raise IOError("No configuration folder detected, please run cblaster config")

    # Check if there's an existing config.ini
    config_ini = config_dir / "config.ini"
    if not config_ini.exists():
        raise IOError("config.ini not found, please run cblaster config")

    # Read in the configuration values and return the parser
    config_parser = configparser.ConfigParser()
    config_parser.read(config_ini)
    return config_parser


def write_config_file(**kwargs):
    LOG.info("Starting config module")

    if not kwargs:
        raise IOError("No arguments given")

    config_dir = get_config_dir()

    # Make directory, if it doesn't exist
    if not config_dir.is_dir():
        LOG.info("No config directory found, creating %s", config_dir)
        config_dir.mkdir(parents=True)

    # Initialise the ConfigParser object
    config_ini = config_dir / "config.ini"
    config_parser = configparser.ConfigParser()

    # Read in any config that already exists
    if config_ini.exists():
        LOG.info("Found existing configuration file, reading")
        config_parser.read(config_ini)
    else:
        LOG.info("No configuration file found, creating")
        config_parser["cblaster"] = {}

    # Set new values
    for key, value in kwargs.items():
        if not value:
            continue
        config_parser["cblaster"][key] = str(value)

    # Write new config file
    LOG.info("Writing configuration to %s", config_ini)
    with config_ini.open("w") as cfg:
        config_parser.write(cfg)

    LOG.info("Done!")
