#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import os
from setuptools import setup, find_packages

# README
with open("README.md", "r") as readme_file:
    readme = readme_file.read()

# REQUIREMENTS
with open("requirements.lock", "r", encoding="utf-8") as fh:
   requirements = [d.strip() for d in fh.readlines() if "#" not in d]

# SCRIPTS
scriptsPath = "gammapysimulator/scripts/"
scripts = [scriptsPath+file for file in os.listdir(scriptsPath)]
for script in scripts:
     print(script)

# SETUP
setup(
    name="gammapysimulator",
    version="0.0.1",
    author="Gabriele Panebianco",
    author_email="gabriele.panebianco@inaf.it",
    description="A package to convert your Jupyter Notebook",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/gabriele-panebianco-inaf/gammapysimulator",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    scripts=scripts,
    package_dir={'gammapysimulator': 'gammapysimulator'},
    include_package_data=True
)
