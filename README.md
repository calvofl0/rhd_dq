rhd_dq
======
*rhd_dq* is a library for reading CO5BOLD UIO files. It includes a python wrapper, *pybold*, and some additional tools. Only the python wrapper is documented hereafter.

Licensing information
---------------------------
Copyright (C) 2019  Flavio Calvo <flavio.calvo@irsol.ch>.

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT license as published by the Open Source
Initiative.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the MIT license for more details.

You should have received a copy of the MIT license along with this program. If
not, see <https://opensource.org/licenses/MIT>.

**Whatever work you do and distribute based on the present code IS HEAVILY
ENCOURAGED to be licensed with a free and open-source software license.**

The files under the 3dparty directory are external 3rd-party software and are
not covered by the present copyright, as well as the files that explicitly
provide copyright information and license. All files outside of the 3dparty
folder and with no explicit copyright notice are part of the present project
and are threfore covered by the present copyright.

Requirements
------------------
The *pybold* python module was tested with both Python 2.7 and Python 3, and with NumPy 1.11.0 and above. *f2py* is required at compile time, and many distributions provide it as a side package *numpy-develop*.

The sources of the CO5BOLD code are also required.

Compilation
---------------
First identify the path to the CO5BOLD source files on your computer and check
the configuration in <CO5BOLD_SRC>/hd/rhd/<YOUR_MACHINE>/.configurerc. The
CO5BOLD code must compile with no errors on your machine, and the F90_POSTFLAGS
variable in the configuration file must contain the compilation flag
"-Drhd_box_arrays01=0".

You can find additional information about CO5BOLD configuration and compilation
in <CO5BOLD_SRC>/README.

Compile the pybold module with
```shell
./configure ARCH=<YOUR_MACHINE> F90_BASEPATH=<CO5BOLD_SRC> COBOLD_3DPARTY=no
make pybold
```
Compilation went successful if the last output line is

```
...
rm -rf "rhd_dq_module.a"
```

and if *_pybold.so* was created. Make sure that both *pybold.py* and *_pybols.so* are in a place where Python can find them.

The pybold interface
--------------------------

In Python, start importing the *pybold* module and creating a model structure. You can easily get help on the model structure:

```python
import pybold
model = pybold.uio_struct()
help(model)
```

The contents of an UIO file can be loaded into the model structure and displayed with

```python
model.load('myfile.uio')
model
```

For instance, you could load a parameter file and display its contents with

```python
model.load('myfile.uio')
model
```

Note however that if the UIO file contains several atmosphere snapshots, only the last one is loaded. If you need to process all the snapshots of an UIO file, say a MEAN file *rhd.mean*, you might prefer to use

```python
model.open('rhd.mean')
model.header        # load the header
model.next          # load the first snapshot
# ...do your processing...
model.next
# ...do your processing...
# ...
model.close
```

Derived quantities with the pybold interface
------------------------------------------------------

Assuming you have a parameter file *rhd.par* with the corresponding files for the equation of state and for opacities, and you are interested in the second snapshot of the FULL file *rhd.full*, you can load it into the model structure using

```python
model.load('rhd.full', 'rhd.par', 2)
```

You will find a new substructure *dq* with all *d*erived *q*uantities:

```python
model.dq
model.dq.T          # compute temperature
```

Derived quantities are computed (and corresponding memory allocated) on an on-request basis.

It is also possible to use the open statement with parameter files:

```python
model.open('rhd.full', 'rhd.par')
model.header        # load the header into the model structure
model.next          # load the first snapshot
model.next          # load the second snapshot
model.close         # close `rhd.full'
```
