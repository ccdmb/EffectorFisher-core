from setuptools import setup, find_packages

setup(
    name="EffectorFisher-core",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.0.0",
    ],
    author="Kristina K. Gagalova",
    author_email="kristina.gagalova@curtin.edu.au",
    description="EffectorFisher: Association of disease phenotype with pangenome-derived protein-isoform profiles for improved prediction of fungal pathogenicity effectors",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/muhitulh/EffectorFisher-core",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
