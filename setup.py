from setuptools import find_packages, setup

with open("requirements.txt") as f:
    required = f.read().splitlines()


setup(
    name="FANDAS",
    license="Apache License 2.0",
    version="2.3.0",
    author="Siddarth Narasimhan, Rodrigo Honorato",
    description="Fast Analysis of multidimensional NMR DAta Sets",
    author_email="",
    include_package_data=True,
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[],
    python_requires=">=3.6, <4",
    install_requires=required,
    entry_points={
        "console_scripts": [
            "fandas=fandas.cli:maincli",
        ],
    },
)
