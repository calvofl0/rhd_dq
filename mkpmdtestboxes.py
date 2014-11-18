#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pybold
import numpy as np

model=pybold.pmd_struct()
model.endianness=0
model.int_sz=4
model.db_sz=8
model.periodic=np.array([1, 0])
model.domain_sz=np.array([20., 30., 40.])
model.domain_origin=np.array([1., 2., 3.])
model.dimensions=np.array([2, 3, 4])
model.xc1=np.array([6., 16.])
model.xc2=np.array([7., 17., 27.])
model.xc3=np.array([8., 18., 28., 38.])
model.n_radtheta=2
model.n_radphi=4
model.atomic_module='PORTA Continuum'
model.pmd_comment='The typical small box for testing'
model.module_header='Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed non risus. Suspendisse lectus tortor, dignissim sit amet, adipiscing nec, ultricies sed, dolor. Cras elementum ultrices diam. Maecenas ligula massa, varius a, semper congue, euismod non, mi. Proin porttitor, orci nec nonummy molestie, enim est eleifend mi, non fermentum diam nisl sit amet erat. Duis semper. Duis arcu massa, scelerisque vitae, consequat in, pretium a, enim. Pellentesque congue. Ut in risus volutpat libero pharetra tempor. Cras vestibulum bibendum augue. Praesent egestas leo in pede. Praesent blandit odio eu enim. Pellentesque sed dui ut augue blandit sodales. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia Curae; Aliquam nibh. Mauris ac mauris sed pede pellentesque fermentum. Maecenas adipiscing ante non diam sodales hendrerit.'
model.module_hd_sz=len(model.module_header)
model.data001=np.empty((2,3,4), order='F')
for i in range(2):
	for j in range(3):
		for k in range(4):
			model.data001[i, j, k] =i+2.*j+6.*k
model.data002=-model.data001
model.pmd_write('PMD_little-endian_int4_db8.pmd')
model.endianness=1
model.pmd_write('PMD_big-endian_int4_db8.pmd')
model.endianness=0
model.int_sz=8
model.db_sz=16
model.pmd_write('PMD_little-endian_int8_db16.pmd')
