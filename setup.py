from setuptools import setup, find_packages, Extension
setup(
    name='spkmeansmodule',
    version='0.1.0',
    author='Bar Jacob & Avital Gendelev',
    author_email="",
    description="full spectral kmeans algorithm",
    install_requires=['invoke'],
    packages=find_packages(),

    license='GPL-2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: Hebrew',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',

    ],
    ext_modules=[
        Extension(
            'spkmeansmodule',
            ['spkmeansmodule.c', 'spkmeans.c', 'kmeans.c'],
        ),
    ]

)
