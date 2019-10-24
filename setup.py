from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()

setup(
    name="clusterBLASTer",
    author="Cameron Gilchrist",
    version="0.0.0",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gamcil/clusterblaster",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["requests"],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["cblaster=clusterblaster.main:main"]},
)
