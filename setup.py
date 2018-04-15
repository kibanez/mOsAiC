from setuptools import setup, find_packages

setup(
    name='mOsAiC',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/kibanez/mosaic',
    install_requires=[
        'HTSeq==0.9.1',
        'PyVCF==0.6.8',
        'pandas==0.22.0',
        'scipy',
        'numpy==1.14.2'
    ],
    license='',
    classifiers=[
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Bioinformatics Data science',
        'Topic :: Software Development :: Mosaicism Bioinformatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7.10',
        'Programming Language :: Python :: 3.6.5',
    ],
    keywords='sample tracking metagenomics microbiotica',
    author='kristina ibanez garikano',
    author_email='kristina.ibanez-garikano@genomicsengland.co.uk',
    description='Detection of somatic mosaicism variants in NGS'
)