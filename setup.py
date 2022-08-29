import setuptools

setuptools.setup(
    name="brixs",
    version="0.1",
    author="Carlos Galdino, Thiago Mori, Tulio Rocha",
    author_email="galdino@ifi.unicamp.br, thiago.mori@lnls.br, tulio.rocha@lnls.br",
    description="Python package for analysis of XAS/RIXS spectra.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/cwgaldino/brixs",
    packages=['brixs', 'brixs.backpack', 'brixs.file_reading'],
    packages_dir={'brixs': './brixs'},
    py_modules=['brixs'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
