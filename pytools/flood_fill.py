#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

# Flood fill algorithm (selects a connected region)
# Implemented with periodic boundary conditions and with the possibility of
# restraining to pixels behind a mask, also records border
def flood_fill(arr, centre, mask=None):
	# Start from empty queue
        Q = []
        if mask is None: mask=np.ones_like(arr, dtype=bool)
        if not arr[centre]:
                return np.zeros_like(arr, dtype=bool), np.zeros_like(arr, dtype=bool)
        arr0 = np.zeros_like(arr, dtype=bool)
        border = np.zeros_like(arr, dtype=bool)
        s = np.shape(arr)
	# Add starting point to queue
        Q.append(centre)
	# Loop over points in queue
        while len(Q) > 0:
		# Remove from queue if outside of domain or already analysed
                if not arr[Q[0]] or not mask[Q[0]] or arr0[Q[0]]:
                        Q.pop(0)
                        continue
                i_min = Q[0][0]
                i_max = Q[0][0]
		# Get index of furthest pixel to the west in domain
                for i in range(Q[0][0]-1,Q[0][0]-s[0],-1):
                        cell = (i%s[0], Q[0][1])
                        if arr[cell] and mask[cell] and not arr0[cell]:
                                i_min = i
                        else: break
		# Record west border (nw and sw included)
                for cell in [((i_min-1)%s[0], (Q[0][1]-1)%s[1]), ((i_min-1)%s[0], Q[0][1]), ((i_min-1)%s[0], (Q[0][1]+1)%s[1])]:
	                if not arr[cell] or not mask[cell]:
				border[cell] = True
		# Get index of furthest pixel to the east in domain
                for i in range(Q[0][0]+1,Q[0][0]+s[0]):
                        cell = (i%s[0], Q[0][1])
                        if arr[cell] and mask[cell] and not arr0[cell]:
                                i_max = i
                        else: break
		# Record east border (ne and se included)
                for cell in [((i_max+1)%s[0], (Q[0][1]-1)%s[1]), ((i_max+1)%s[0], Q[0][1]), ((i_max+1)%s[0], (Q[0][1]+1)%s[1])]:
	                if not arr[cell] or not mask[cell]:
				border[cell] = True
		# Add to the queue all pixels north and south of the
		# selected segment and mark the segment
                for i in range(i_min, i_max+1):
                        arr0[i%s[0], Q[0][1]] = True
                        north = (i%s[0], (Q[0][1]+1)%s[1])
                        nw = ((i-1)%s[0], (Q[0][1]+1)%s[1])
                        ne = ((i+1)%s[0], (Q[0][1]+1)%s[1])
                        south = (i%s[0], (Q[0][1]-1)%s[1])
                        sw = ((i-1)%s[0], (Q[0][1]-1)%s[1])
                        se = ((i+1)%s[0], (Q[0][1]-1)%s[1])
                        if arr[north] and mask[north] and not arr0[north]:
				Q.append(north)
                        if arr[south] and mask[south] and not arr0[south]:
				Q.append(south)
			# Record north and south borders
			for cell in [north, south, ne, nw, se, sw]:
	                        if not arr[cell] or not mask[cell]:
					border[cell] = True
        return arr0, border
