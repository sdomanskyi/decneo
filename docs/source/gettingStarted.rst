.. _getting-started:

**Getting Started**
===================

These instructions will get you a copy of the project up and running on your machine for data analysis, development or testing purposes.

**Installation**
----------------

Install of the latest release of ``decneo``:

.. code-block:: bash

    $ pip install decneo

For detailed instructions and other ways to install ``decneo`` as well as
list of optional packages and instructions on how to install them see
**Prerequisites** section at https://github.com/sdomanskyi/decneo


**Loading the package**
-----------------------

In your script import the package:

.. code-block:: python

	from decneo.analysisPipeline import Analysis

Create an instance of ``class decneo``. Here, for simplicity, we use Default parameter values:

.. code-block:: python

	an = Analysis()
