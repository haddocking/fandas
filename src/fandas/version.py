from importlib.metadata import version as package_version

VERSION = package_version("fandas")

v_major, v_minor, v_patch = VERSION.split(".")
