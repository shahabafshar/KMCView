#!/usr/bin/env python3
"""
Setup configuration for Zacros Molecular Evolution Viewer
"""

from setuptools import setup, find_packages
import os

# Read long description from README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements from requirements.txt
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

# Filter out development dependencies
install_requires = [req for req in requirements if not any(dev in req for dev in ["pytest", "sphinx", "flake8", "black", "isort", "mypy"])]

setup(
    name="zacros-molecular-viewer",
    version="1.0.0",
    author="Zacros Visualization Team",
    author_email="contact@zacros-viewer.com",
    description="Interactive visualization tools for Zacros kinetic Monte Carlo simulation data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/zacros-molecular-viewer",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/zacros-molecular-viewer/issues",
        "Source": "https://github.com/yourusername/zacros-molecular-viewer",
        "Documentation": "https://github.com/yourusername/zacros-molecular-viewer/wiki",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=install_requires,
    extras_require={
        "dev": ["pytest>=6.0.0", "pytest-cov>=2.12.0", "flake8>=3.9.0", "black>=21.0.0", "isort>=5.9.0", "mypy>=0.910"],
        "docs": ["sphinx>=4.0.0", "sphinx-rtd-theme>=1.0.0"],
    },
    entry_points={
        "console_scripts": [
            "zacros-matplotlib-viewer=fixed_matplotlib_viewer:main",
            "zacros-web-viewer=fixed_web_viewer:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.dat"],
    },
    zip_safe=False,
    keywords="zacros, molecular dynamics, visualization, monte carlo, simulation, chemistry, physics",
) 