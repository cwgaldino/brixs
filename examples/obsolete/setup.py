import setuptools

setuptools.setup(
    name="brixs",
    version="0.9.6",
    author="Carlos Galdino, Thiago Mori, Felipe Custodio, Tulio Rocha",
    author_email="carlos.galdino@psi.ch, thiago.mori@lnls.br, felipe.custodio@lnls.br, tulio.rocha@lnls.br",
    description="python package for processing and analysis of XAS and RIXS spectra ",
    long_description=open('README.rst').read(),
    long_description_content_type="text/x-rst",
    url="https://github.com/cwgaldino/brixs",
    packages=['brixs', 'backpack', 'beamlines'],
    packages_dir={'brixs': './brixs', 'backpack': './backpack', 'beamlines': './beamlines'},
    py_modules=['brixs'],
    install_requires=['matplotlib>=3.4', 'numpy>=1.15'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
