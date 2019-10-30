import re
from pathlib import Path

from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()


def get_version():
    """Get version number from __init__.py"""
    version_file = Path("clusterblaster/__init__.py").read_text()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


setup(
    name="clusterBLASTer",
    author="Cameron Gilchrist",
    version=get_version(),
    description="",
    long_description=long_description,
    license="MIT",
    url="https://github.com/gamcil/clusterblaster",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: MIT",
    ],
    install_requires=["requests"],
    tests_require=["pytest-cov", "pytest-mock", "requests-mock"],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["cblaster=clusterblaster.main:main"]},
)
