from setuptools import setup, find_packages

setup(
    name="FANDAS",
    license="Apache License 2.0",
    version="2.0.1",
    author="Siddarth Narasimhan",
    description="Fast Analysis of multidimensional NMR DAta Sets",
    author_email="",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[],
    python_requires=">=3.6, <4",
    install_requires=["numpy"],
    entry_points={
        "console_scripts": [
            "fandas=fandas.cli:maincli",
        ],
    },
)
