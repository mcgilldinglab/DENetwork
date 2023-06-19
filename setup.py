from setuptools import setup

setup(  name='DENetwork',
        version='1.0',
        description='Unveiling Regulatory and Signaling Networks Behind Differentially Expressed Genes',
        author='Ting-Yi Su',
        author_email='ting-yi.su@mail.mcgill.ca',
        url="https://github.com/mcgilldinglab/DENetwork",
        license='MIT',
        packages=['DENetwork'],
        install_requires=['numpy>=1.19.5','networkx>=2.5.1','pandas>=1.1.5','matplotlib>=3.3.4','scipy>=1.5.4','kneed>=0.7.0','statsmodels>=0.12.2','openpyxl>=3.1.2','rpy2>=3.4.4'],
        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
        ],
        )