from setuptools import Extension, setup

module = Extension("symnmf", 
                   sources=['symnmfmodule.c', 'utils.c', 'sym.c', 'diagonal.c', 'norm.c', 'symnmf.c'],
                   extra_compile_args=['-g'] 
)
setup(name='symnmf',
     version='1.0',
     description='Python wrapper for SymNMF in C',
     ext_modules=[module]
    )

