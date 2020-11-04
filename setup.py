import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="brixs",
    version="0.1",
    author="Carlos Galdino",
    author_email="galdino@ifi.unicamp.br",
    description="Python package for analysis of RIXS spectra.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cwgaldino/brixs",
    packages=['brixs'],
    packages_dir={'brixs': './brixs'},
    py_modules=['brixs',
                'brixs.arraymanip',
                'brixs.figmanip',
                'brixs.filemanip',
                'brixs.intermanip',
                'brixs.model_functions',
                'brixs.simulate'],
    # package_data={'brixs': ['./default_parameters.dat']
    #                 },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
