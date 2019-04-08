#!/usr/bin/env python

from distutils.core import setup

setup(name='REDItools',
	version='1.0.4',
	description='Python Scripts for RNA editing detection by RNA-Seq data',
	author='Ernesto Picardi',
	author_email='ernesto.picardi@gmail.com',
	url='http://code.google.com/p/reditools/',
	scripts=['reditools/REDItoolBlatCorrection.py',
	'reditools/REDItoolDenovo.py',
	'reditools/REDItoolDnaRna.py',
	'reditools/REDItoolKnown.py',
	'reditools/AnnotateTable.py',
	'reditools/FilterTable.py',
	'reditools/SearchInTable.py',
	'reditools/selectPositions.py',
	'reditools/GFFtoTabix.py',
	'reditools/SortGFF.py',
	'reditools/SortTable.py',
	'reditools/TableToGFF.py',
	'reditools/tableToTabix.py',
	],
	license='LICENSE.txt',
	classifiers=[
          'Intended Audience :: Computational biologists',
          'License :: OSI Approved :: MIT',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          ],
	long_description=open('README').read(),
	platforms=['Linux','Unix','MacOS']
)

