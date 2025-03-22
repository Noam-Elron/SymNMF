from setuptools import Extension, setup

module = Extension("symnmf", sources=['symnmfmodule.c'])
setup(name='symnmf',
     version='1.0',
     description='Python wrapper for SymNMF in C',
     ext_modules=[module]
    )

