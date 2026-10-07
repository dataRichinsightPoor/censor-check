# Bundled data

`chembl-six-targets-journal-ic50.csv` contains 31,236 IC50 records retrieved from the ChEMBL web services (https://www.ebi.ac.uk/chembl/) on October 7, 2026. The records cover six single-protein targets: VEGFR2 (CHEMBL279), ABL1 (CHEMBL1862), AKT1 (CHEMBL4282), BACE1 (CHEMBL4822), thrombin (CHEMBL204) and HDAC1 (CHEMBL325).

They are restricted to records whose source is the scientific literature (ChEMBL source identifier 1) with standard units of nM and a positive standard value. The database release number could not be confirmed, because the service's status endpoint returned an error during retrieval.

Columns: `document` (ChEMBL document ID), `molecule` (ChEMBL molecule ID), `target`, `relation` (standard relation), `value_nM` (standard value), `units`, `validity` (data validity comment), `year`.

ChEMBL data are made available under the Creative Commons Attribution-ShareAlike 3.0 Unported license (https://creativecommons.org/licenses/by-sa/3.0/). This derived subset is distributed under the same license. Please cite ChEMBL: Zdrazil B, et al. Nucleic Acids Research. 2024;52(D1):D1180–D1192. https://doi.org/10.1093/nar/gkad1004
