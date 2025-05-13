from setuptools import setup, find_packages

with open("README.md", "r") as file :
    doc = file.read()

setup(
    name="md_manager",
    version="2.3.1",
    description="A python module that allow molecular dynamics data analysis based on pandas DataFrames.",
    long_description=doc,
    long_description_content_type="text/markdown",

    packages=find_packages(),
    install_requires=["mdanalysis>=2.8"],
    extras_requires={"dev":["twine>=4.0.2"]},
    python_requires=">=3.8",

    url="https://github.com/NicolasPetiot/md_manager",

    author="NicolasPetiot",
    author_email="nicolaspetiot2710@hotmail.fr",

    license="AGPL",
    classifiers=[
        "Programming Language :: Python :: 3.8"
    ]
)

