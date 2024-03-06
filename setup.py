from setuptools import setup, find_packages

with open("README.md", "r") as file :
    doc = file.read()

setup(
    name="md_manager",
    version="0.0.1",
    description="A python module that allow molecular dynamics data analysis based on pandas DataFrames.",
    long_description=doc,

    packages=find_packages(),
    install_requires=["pandas", "numpy", "numba", "scipy", "networkx"],

    url="https://github.com/NicolasPetiot/md_manager",

    author="NicolasPetiot",
    author_email="nicolaspetiot2710@hotmail.fr",
)
