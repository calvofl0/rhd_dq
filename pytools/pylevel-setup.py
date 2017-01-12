#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
		name         = 'pylevel',
		version      = '0.1.0',
		description  = 'This module provides the 3D level function',
		author       = 'Flavio Calvo',
		author_email = 'flavio.calvo@irsol.ch',
		ext_modules  = [Extension("pylevel", ["pylevel.c"], libraries=['m'], extra_compile_args=[""])],
		include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
)
