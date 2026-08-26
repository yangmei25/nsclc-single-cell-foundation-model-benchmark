# Third-party data, software, and publications

The repository-level MIT License applies only to the original code and
documentation authored for this demo. It does not relicense third-party data,
publications, software, pretrained models, or model weights. Those materials
remain subject to their respective terms.

## GSE205335 source data

This project reanalyzes selected, publicly available processed single-cell
RNA-sequencing data and metadata from NCBI Gene Expression Omnibus accession
[GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335):

> Single-cell transcriptome profiles of tumor tissues from lung cancer
> patients receiving immune checkpoint inhibitors.

NCBI places no restrictions on the use or distribution of GEO data, while
noting that submitters may retain or assert rights in submitted material. See
the [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
Users of this repository remain responsible for following the terms attached
to the source record and for citing the accession and associated publication.

Raw or processed source files downloaded from GEO are intentionally excluded
from this repository. Small cohort tables and figures committed here are
derived analysis artifacts; they are not substitutes for the authoritative
GEO record.

## Associated publication

Kim N, Park S, Jo A, et al. (2024). *Unveiling the influence of tumor and
immune signatures on immune checkpoint therapy in advanced lung cancer*.
eLife 13:RP98366. <https://doi.org/10.7554/eLife.98366.3>

Copyright © 2024 Kim, Park et al. The article is distributed under the
[Creative Commons Attribution License](https://creativecommons.org/licenses/by/4.0/),
which permits reuse with attribution. This repository does not redistribute
the article PDF or copy its published figures. Its figures and conclusions
were generated independently from the selected data subset, and changes in
cohort selection and analysis are identified in the README and notebooks.

## scGPT

This project uses the external [scGPT](https://github.com/bowang-lab/scGPT)
Python package and a pretrained scGPT checkpoint to compute frozen cell
embeddings. The scGPT source repository is distributed under the
[MIT License](https://github.com/bowang-lab/scGPT/blob/main/LICENSE).

Cui H, Wang C, Maan H, et al. (2024). *scGPT: toward building a foundation
model for single-cell multi-omics using generative AI*. Nature Methods 21,
1470–1480. <https://doi.org/10.1038/s41592-024-02201-0>

scGPT source code and pretrained model weights are not redistributed by this
repository. Users must obtain them from their official source and comply with
any terms supplied with the particular checkpoint. Merely mentioning or
depending on scGPT does not place scGPT under this repository's license.

## Other dependencies

Scanpy, AnnData, PyTorch, pandas, NumPy, scikit-learn, and the other packages
listed in the environment files are third-party projects. Installation through
the environment specifications is subject to each package's own license.

No affiliation with or endorsement by the source-study authors, NCBI, eLife,
or the scGPT developers is claimed.
