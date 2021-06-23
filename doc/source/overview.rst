Overview
========

Notations
---------
* :math:`S^{\star}_t` : segmentation at time :math:`t`

* :math:`S^{\star}_{t+\delta t \leftarrow t} = S^{\star}_t \circ \mathcal{T}_{t \leftarrow t+\delta t}` : segmentation à :math:`t` projetée dans :math:`I_{t+\delta t}`. C'est une nouvelle notation par rapport à :cite:`guignard:tel-01278725`, mais l'image est pas mal utilisée.

*  :math:`\tilde{S}_{t+1}` : segmentation de :math:`I_{t + \delta t}` par ligne de partage des eaux avec les graines :math:`S^e_{t+1 \leftarrow t} = S^e_t \circ \mathcal{T}_{t \leftarrow t+\delta t}`

* :math:`\hat{S}_{t+1}` : segmentation de :math:`I_{t + \delta t}` par ligne de partage des eaux avec les graines :math:`\mathrm{Seeds}_{t+1}`. Cette image peut être amenée à changer

* :math:`S^e_t` : cellules de :math:`S^{\star}_t` érodées (pour servir de graines)

* :math:`S^e_{t+\delta t \leftarrow t} = S^e_t \circ \mathcal{T}_{t \leftarrow t+\delta t}` : cellules de :math:`S^{\star}_t` érodées (pour servir de graines) puis projetées dans `I_{t+\delta t}`.

* :math:`\mathrm{Seeds}_{t+1}` : image des graines calculées avec les paramètres optimaux. Cette image peut être amenée à changer

* :math:`\mathcal{T}_{t \leftarrow t+\delta t}` : transformation non-linéaire permettant de rééchantillonner :math:`I_{t}` sur  :math:`I_{t + \delta t}`

