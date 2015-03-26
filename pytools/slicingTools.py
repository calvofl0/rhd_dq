#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def vorticity(v1, v2, delta=1.):
	n1, n2 = np.shape(v1)
	assert(np.shape(v1) == np.shape(v2))
	dv1dx2 = (v1[:,range(1,n1)+[0]]-v1)/delta
	dv2dx1 = (v2[range(1,n2)+[0],:]-v2)/delta
	return dv1dx2 - dv2dx1
