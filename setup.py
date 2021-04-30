import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="LyaPlotter", # Replace with your own username
    version="0.0.1",
    author="César Ramírez",
    author_email="cramirez@ifae.es",
    description="Lyman Alpha mock pipeline sims tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cramirezpe/LyaPlotter/",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'LyaPlotter_get_pixel_neighbours = LyaPlotter.scripts.get_pixel_neighbours:main',
            'LyaPlotter_master_to_qso_cat = LyaPlotter.scripts.tools_scripts:master_to_qso_cat_script',
            'LyaPlotter_test_debug = LyaPlotter.tools:test_debug_mode',
            'LyaPlotter_colore_to_drq = LyaPlotter.scripts.tools_scripts:colore_to_drq_script',
            'LyaPlotter_trim_catalog_into_pixels = LyaPlotter.scripts.tools_scripts:trim_catalog_into_pixels_script',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)   