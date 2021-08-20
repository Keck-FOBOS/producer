
.. include:: include/links.rst

Installation
============

Clone the repo
--------------

To download the software and associated data, clone the GitHub repo by
executing:

    .. code-block:: console

        git clone https://github.com/Keck-FOBOS/producer.git

This will create a ``producer`` directory in the current directory.

Install Python 3
----------------

The FOBOS Producer is supported for Python 3 only. To install Python,
you can do so along with a full package manager, like `Anaconda`_, or
you can install python 3 directly from `python.org`_.


Install from source
-------------------

The preferred method to install ``producer`` and ensure its dependencies are
met is to, from the top-level directory, run:

.. code-block:: console

    pip install -e .

This approach is preferred because it eases uninstalling the code:

.. code-block:: console
    
    pip uninstall fobos-producer

Installation in this way should also mean that changes made to the code
should take immediate effect when you restart the calling python
session.

----

Test your installation
----------------------

To test the installation, use ``pytest``:

.. code-block:: console

    cd producer/tests
    pytest . -W ignore

Problems?
---------

We have limited support to offer installation help. However, if you
have problems, particularly those that you think may be a more
general problem, please `Submit an issue`_.

