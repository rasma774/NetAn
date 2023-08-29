NetAn
================

About NetAn
-------------
Gene annotation enrichment analysis is the gold standard to study the biological context of a set of genes, but available tools often overlooks important network properties of the gene regulatory system. We introduce the Network Annotation Enrichment package, NetAn, which takes a list of genes and uses network-based approaches such as network clustering and inference of closely related genes to include local neighbours. Using NetAn, we exemplify how these approaches can be used to better identify relevant annotations in human gene sets. 

Examples of how to use NetAn are found below, and in a Jupyter notebook under ./examples/

Installation
============
The package can be used under the GNU GENERAL PUBLIC LICENSE V3

To install 
```consol
pip install .
```

To test the installation, use e.g. pytest 

```consol
pip install pytest
pytest test/test.py 
```


Python Version Supported/Tested
-------------------------------
- Python 3.11


Usage:
======
Here are the basic usages of the TFTenricher. A more comprehensive example can e found as a Jupyter notebook in the ./doc/ folder

In python:
```python
from NetAn import NetAn

test_genes = ['Gata3', 'STAT2', 'INSR']
enr = NetAn(test_genes)

# To get the results of a clustering run
# The gene_set parameter sets wheather  
cluster_res = enr.enrich(gene_set='start_genes',
                         cluster=True,
                         n_clusters=2,
                         )


# to print the inferred target genes
print(enr.target_genes)
```

Contributor:
=============

Rasmus Magnusson: Development of the package.

Current Members in the Project
------------------------------
- @rasma774

References & how to cite
======================
n/a
