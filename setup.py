import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setuptools.setup(
    name="SEDFit",
    version="1.0.0",
    author="Marina Kounkel",
    author_email="marina.kounkel@unf.edu",
    description="Performs a model SED fitting to the stellar fluxes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mkounkel/SEDFit",
    packages=['SEDFit'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['astropy','astroquery','numpy','scipy','dust_extinction','dustmaps','matplotlib','tqdm','tensorflow>2.16'],
    package_data={
        'SEDFit': ['sed_coelho.p','sed_kurucz.p','sed_phoenix.p','sed_btsettl.p','*.dat','quality.p']
    },
    python_requires='>=3.6',
)
