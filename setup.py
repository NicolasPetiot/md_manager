from setuptools import setup, find_packages

with open("README.md", "r") as file :
    doc = file.read()

setup(
    name="md_manager",
    version="0.0.1",
    description="",
    package_dir={"":"app"},
    packages=find_packages(where="app"),
    long_description=doc,
    #long_description_type="text/markdown",
    url="https://github.com/NicolasPetiot/md_manager",
    author="NicolasPetiot",
    author_email="nicolaspetiot2710@hotmail.fr",
    #license="",
    #classifiers=[],
    install_requires=[
        "pandas >= 2.1.1",
        "numpy >= 1.26.2",
        "numba >= 0.58.1",
        "scipy >= 1.11.4",
        "networkx >= 3.1"
    ],
    extras_require={
        "dev": ["twine>=4.0.2"]
    },
    python_requires="<3.12"

)

