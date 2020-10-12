**Getting Started**
===================

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

**Installation**
----------------

Install of the latest release of ``scRegulation``:

.. code-block:: bash

    $ pip install scRegulation

For detailed instructions and other ways to install ``scRegulation`` as well as
list of optional packages and instructions on how to install them see
**Prerequisites** section at https://github.com/sdomanskyi/scRegulation


**Loading the package**
-----------------------

In your script import the package:

.. code-block:: python

	from scRegulation.analysisPipeline import Analysis

Create an instance of ``class scRegulation``. Here, for simplicity, we use Default parameter values:

.. code-block:: python

	an = Analysis()
