from setuptools import setup, find_packages

setup(
    name="cosmo_lensing",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.24",
        "scipy>=1.10",
        "astropy>=5.0",
        "matplotlib>=3.6",
        "treecorr>=4.3",
        "healpy>=1.16",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=23.0",
            "isort>=5.12",
        ],
        "workflow": [
            "snakemake>=7.0",
        ],
    },
    python_requires=">=3.9",
    author="Laurent Magri-Stella",
    description="Weak gravitational lensing analysis pipeline for cosmological N-body simulations",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3.9",
    ],
)
