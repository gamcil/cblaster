import re
from pathlib import Path

from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()


def get_version():
    """Get version number from __init__.py"""
    version_file = Path(__file__).resolve().parent / "cblaster" / "__init__.py"
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read_text(), re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


setup(
    name="cblaster",
    author="Cameron Gilchrist",
    version=get_version(),
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/gamcil/cblaster",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
        "biopython",
        "requests",
        "numpy",
        "scipy",
        "PySimpleGUI",
        "Biopython",
        "clinker>=0.0.15",
        "gffutils",
    ],
    tests_require=["pytest", "pytest-cov", "pytest-mock", "requests-mock"],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["cblaster=cblaster.main:main"]},
    package_data={"cblaster": ["plot/*"]},
)
