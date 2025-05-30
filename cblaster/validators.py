# Validator functions for argument parsers

import argparse
import os

from multiprocessing import cpu_count
from pathlib import Path


def full_path(file_path, *acces_modes, dir=False):
    """Test if a file path or directory exists and has the correct permissions and create a full path

    For reading acces the file has to be pressent and there has to be read acces. For writing acces the directory with
    the file has to be present and there has to be write acces in that directory.

    Args:
        file_path (str): relative or absoluete path to a file
        acces_modes (List): a list of integers of acces modes for which at least one should be allowed
        dir (bool): if the path is to a directory or not
    Returns:
        A string that is the full path to the provided file_path
    Raises:
        argparse.ArgumentTypeError when the provided path does not exist or the file does not have the correct
        permissions to be accessed
    """
    full_file_path = Path(file_path).absolute().resolve()
    failed_path = failed_acces = False
    for acces_mode in acces_modes:
        if not dir and full_file_path.is_file() or (acces_mode == os.W_OK and full_file_path.parent.is_dir()):
            if os.access(full_file_path, acces_mode) or \
                    (acces_mode == os.W_OK and os.access(full_file_path.parent, acces_mode)):
                return str(full_file_path)
            else:
                failed_acces = True
        elif dir and full_file_path.is_dir():
            if os.access(full_file_path, acces_mode):
                return str(full_file_path)
        else:
            failed_path = True
    if failed_path:
        raise argparse.ArgumentTypeError(f"Invalid path: '{file_path}'.")
    elif failed_acces:
        raise argparse.ArgumentTypeError(f"Invalid acces for path: {file_path}.")


def full_database_path(database, *acces_modes):
    """Make sure the database path is also correct, but do not check when providing one of the NCBI databases

    Args:
        database (str): a string that is the path to the database creation files or a NCBI database identifier
        acces_modes (List): a list of integers of acces modes for which at least one should be allowed
    Returns:
        a string that is the full path to the database file or a NCBI database identifier
    """
    if database not in NCBI_DATABASES:
        try:
            return full_path(database, *acces_modes)
        except argparse.ArgumentTypeError as e:
            raise type(e)(str(e) + f" Or use one of the following databases {', '.join(NCBI_DATABASES)}"
                                   f" when running in remote mode")
    return database


def max_cpus(value):
    """
    Ensure that the cpu's do not go above the available amount. Setting to high cpu's will crash database creation
    badly

    Args:
        value(int): number of cpu's as provided by the user
    Returns:
        value as an integer with 1 <= value <= multiprocessing.cpu_count()
    """
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid 'int' value: f")
    value = max(1, min(value, cpu_count()))
    return value


def int_or_pair(astring):
    """Parses arguments that can be int or (int, int). """
    if astring.isdigit():
        return int(astring)
    if "," not in astring:
        raise argparse.ArgumentTypeError(f"{astring} is not an number, or comma-separated list of numbers")
    try:
        left, right = astring.split(",")
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{astring}' invalid, more than two values given")
    try:
        left = int(left)
        right = int(right)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{astring}' invalid, expected two numbers")
    return (left, right)
