from setuptools import setup, find_packages
from MODtools.version import version

setup(
    name='modtools',
    version=version(),
    packages=find_packages(exclude=['CGRtools', 'CGRtools.cli', 'CGRtools.files']),
    url='https://github.com/stsouko/MODtools',
    license='AGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    description='Modeler tools',
    scripts=['modelbuilder.py'],
    package_data={'': ['unwanted.elem', 'standardrules_dragos.rules']},
    requires=['CGRtools', 'networkx', 'periodictable', 'pandas', 'dill', 'sortedcontainers', 'sklearn', 'numpy',
              'requests'],
    long_description='Modeler tools distributive. include ISIDA Fragmentor and EED python wrappers',

    keywords="tools modeler cli ISIDA Framentor EED SVM IAP",
    classifiers=['Environment :: Console',
                 'Intended Audience :: End Users/Desktop',
                 'Intended Audience :: Developers',
                 ('License :: OSI Approved :: GNU Affero General Public License'
                  ' v3 or later (AGPLv3+)'),
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.2',
                 'Programming Language :: Python :: 3.3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 ]
)
