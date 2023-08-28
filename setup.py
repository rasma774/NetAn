import setuptools

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2022 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'



with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NetAn",
    version="0.01",
    author="Rasmus Magnusson",
    author_email="rasma774@gmail.com",
    description="A package for annotation enrichment analysis from a network perspective",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE V3",
        ],
    python_requires='>=3.11.3',
    install_requires=[
        'scikit-learn>=1.1.0',
        'statsmodels>=0.13.5',
        'matplotlib>=3.7.1',
        'numpy>=1.24.2',
        'scipy>=1.10.1',
        'pandas>=2.0.0',
        'networkx>=3.1',
        ]
)
