import setuptools

setuptools.setup(
    name="brixs",
    version="0.1",
    author="Carlos Galdino, Thiago Mori, Tulio Rocha",
    author_email="galdino@ifi.unicamp.br, thiago.mori@lnls.br, tulio.rocha@lnls.br",
    description="Python package for analysis of RIXS spectra.",
    long_description=open('README.rst').read(),
    long_description_content_type="text/x-rst",
    url="https://github.com/cwgaldino/brixs",
    packages=['brixs'],
    packages_dir={'brixs': './brixs'},
    py_modules=['brixs'],
    # package_data={'brixs': ['./default_parameters.dat']
    #                 },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
