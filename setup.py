from distutils.core import setup

setup(
    name='Seriation',
    version='v2.0',
    packages=['app'],
    url='http://www.lipolab.org/seriation',
    license='Apache Software License, Version 2.0',
    author='Carl P. Lipo',
    author_email='carl.lipo@csulb.edu',
    description='A package for constructing iteratively deterministic frequency seriation solutions.',
        # Include additional files into the package
    include_package_data=True,
        # Dependent packages (distributions)
    install_requires=[
        "flask",
    ],
)

