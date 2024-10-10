from setuptools import setup, find_packages

setup(
    name="pynlcd",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "gdal",
        "numpy",
    ],
    entry_points={
        "console_scripts": [
            "pynlcd=pynlcd.nlcd:main",
        ],
    },
)