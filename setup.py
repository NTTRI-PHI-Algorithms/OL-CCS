from setuptools import setup, find_packages

setup(
    name="OL-CCS",
    version="1.2",
    description="Open-Loop Coherent Compressed Sensor",
    author="M. D. Sudeera H. Gunathilaka",
    author_email="mastiyage.s.aa@m.titech.ac.jp",
    packages=find_packages(),
    install_requires=[
        "numpy==1.26.0",
        "pandas==2.0.3",
        "matplotlib==3.7.2",
        "torch==2.0.1",
        ],
)