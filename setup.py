#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2023/10/7 16:09
# @Author  : lvtongxuan
# @File    : setup.py
# @Software: PyCharm
# @Email   : lvtongxuan@genomics.cn
import setuptools
from wheel.bdist_wheel import bdist_wheel

__version__ = "1.0.3"


class BDistWheel(bdist_wheel):
    def get_tag(self):
        return (self.python_tag, "none", "any")


cmdclass = {
    "bdist_wheel": BDistWheel,
}

requirements = open("requirements.txt").readline()

setuptools.setup(
    name="EAGS",
    version=__version__,
    author="lvtongxuan",
    author_email="lvtongxuan@genomics.cn",
    url="https://github.com/STOmics/EAGS.git",
    description="EAGS: efficient and adaptive Gaussian smoothing applied to high-resolved spatial transcriptomics",
    python_requires=">=3.8",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    cmdclass=cmdclass,
)
