import setuptools
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

scripts = glob.glob('bin/*')

setuptools.setup(
    name="LyaPlotter_cramirezpe", # Replace with your own username
    version="0.0.1",
    author="César Ramírez",
    author_email="cramirez@ifae.es",
    description="Lyman Alpha mock pipeline sims tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cramirezpe/LyaPlotter/",
    packages=setuptools.find_packages(),
    scripts = scripts,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)   