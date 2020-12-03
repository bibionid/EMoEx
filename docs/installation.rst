Installation
============

EMotEP source code: `Github <https://github.com/bibionid/EMotEP>`_

To install the pipeline follow these steps:

1. Clone workflow into working directory

.. code-block:: bash

  git clone https://github.com/bibionid/EMotEP.git path/to/workdir
  cd path/to/workdir

2. Edit config and workflow as needed

.. code-block:: bash

  vim config/emotep_config.yaml

3. Execute workflow, deploy software dependencies via conda

.. code-block:: bash

  snakemake -n --use-conda

Please report any prombels with installation via the `issue tracker  <https://github.com/bibionid/EMotEP/issues>`_
