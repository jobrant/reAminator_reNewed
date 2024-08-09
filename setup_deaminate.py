from setuptools import setup, Extension

module = Extension('deaminate',
                   sources = ['deaminate.c'])

setup(name = 'Deaminate',
      version = '1.0',
      description = 'This is a deamination package',
      ext_modules = [module])