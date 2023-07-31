import setuptools

setuptools.setup(
    name="brixs",
    version="0.9",
    author="Carlos Galdino, Thiago Mori, Tulio Rocha, Felipe Custodio",
    author_email="galdino@ifi.unicamp.br, thiago.mori@lnls.br, tulio.rocha@lnls.br, felipe.custodio@lnls.br",
    description="python package for processing and analysis of XAS and RIXS spectra ",
    long_description=open('README.rst').read(),
    long_description_content_type="text/x-rst",
    url="https://github.com/cwgaldino/brixs",
    packages=['brixs'],
    packages_dir={'brixs': './brixs'},
    py_modules=['brixs'],
    install_requires=['matplotlib>=3.4', 'numpy>=1.15', 'scipy>=1.7', 'lmfit>=1.2.2'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
