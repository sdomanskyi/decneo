About
-----

In silico detection of transcriptional regulation genes from single cell transcriptomics


For development and use:
------------------------

> Note: this is not intended for public use yet. Publication is being prepared.


To install scRegulation as a package:

	$ python setup.py install


For development and testing of the documentation locally (on the development machine) install Sphinx by:

	$ pip install -U sphinx

To compile html version of the documentation go to docs/ directory and run:

	$ sphinx-build -E -a -b html ./source ./build

We are utilizing a 3rd party Sphinx extension sphinxcontrib-images extension, allowing to display documentation images in a organized way.