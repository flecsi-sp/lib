.. |br| raw:: html

   <br />

.. _build:

Prerequisites
*************

* **Doxygen** |br|
  Doxygen is used to build the interface documentation. Doxygen is only
  required if you intend to build the documetnation.

* **Sphinx** |br|
  Sphinx is used to build the main web-based documentation. We also
  require the Sphinx RTD Theme. These can be installed on most Linux
  systems using pip. Sphinx is only required if you intend to build the
  documentation.

Building the Documentation (Developers)
+++++++++++++++++++++++++++++++++++++++

The Burton Specialization uses Doxygen for its API reference and Sphinx
for user and developer documentation.

Doxygen can be installed with most Linux package managers or by using
``spack``.  To install Sphinx, you can install ``pip3`` and use it to
install ``Sphinx``, ``recommonmark``, ``sphinx_rtd_theme``, and
``sphinxcontrib-tikz``. Your package manager should also have ``pip3``;
e.g., on Ubuntu, you can install all of these requirements like:

.. code-block:: console

  $ sudo apt install doxygen
  $ sudo apt install python3-pip
  $ pip3 install --user Sphinx
  $ pip3 install --user recommonmark
  $ pip3 install --user sphinx_rtd_theme
  $ pip3 install --user sphinxcontrib-tikz

To enable documentation, do this:

.. code-block:: console

  $ cmake -DENABLE_DOCUMENTATION=ON ..

By default, this will enable Doxygen and Sphinx. Once you have properly
configured FleCSI, you can build the documentation like:

.. code-block:: console

  $ make doc

This will build the main documentation target in your build directory at
``doc/index.html``. You can open this in your browser with
``file:///path/to/your/build/directory/doc/index.html``.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
