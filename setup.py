# Author: Kobie Kirven
# Penn State University
# Date: 3-7-2022

import setuptools
import sys
from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


CURRENT_DIR = Path(__file__).parent
sys.path.insert(0, str(CURRENT_DIR))


setuptools.setup(
    name="minetwister",
    version="0.01",
    packages=["bin"],
    include_package_data=True,
    entry_points={"console_scripts": ["minetwister=bin.__main__:minetwister",],},
    author="Kobie Kirven",
    description="mineTwister: a tool for mining genomes for twister ribozymes",
    install_requires=["setuptools", "biopython", "pandas"],
    python_requires=">=3.5",
)