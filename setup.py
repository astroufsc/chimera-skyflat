from distutils.core import setup

setup(
    name='chimera_skyflat',
    version='0.0.1',
    packages=['chimera_skyflat', 'chimera_skyflat.controllers', 'chimera_skyflat.interfaces'],
    scripts=[],
    url='http://github.com/astroufsc/chimera-skyflat',
    license='GPL v2',
    author='Antonio Kanaan',
    author_email='kanaan@astro.ufsc.br',
    description='Chimera plugin for automated skyflats'
)
