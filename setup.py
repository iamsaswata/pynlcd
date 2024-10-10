from setuptools import setup, find_packages

setup(
    name="pynlcd",
    version="0.1.3",
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


# from setuptools import setup

# setup(
#   name="pynlcd",
#   version="0.1.2",
#   packages=["pynlcd"],
#   install_requires=[
#       "gdal",
#       "numpy",
#   ],
#   entry_points={
#       "console_scripts": [
#           "pynlcd=pynlcd:main",
#       ],
#   },
# )


# from setuptools import setup

# setup(
#   name="pynlcd",
#   version="0.1.2",
#   packages=["pynlcd"],
#   install_requires=[
#       "gdal",
#       "numpy",
#   ],
#   entry_points={
#       "console_scripts": [
#           "pynlcd=pynlcd.nlcd:main",
#       ],
#   },
# )