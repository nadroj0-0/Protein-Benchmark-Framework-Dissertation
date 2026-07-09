# CAFA3 Historical Validation Report

Final status: **PASS**

Status thresholds:
- PASS: every CSV exists, minimum label F1/protein Jaccard/GO-term Jaccard >= 0.9999.
- PASS WITH MINOR DIFFERENCES: every CSV exists and minimum label F1 >= 0.95.
- FAIL: missing required CSVs or larger discrepancies.

## CSV Summary

|file|generated_rows|reference_rows|row_diff|protein_jaccard|go_term_jaccard|matrix_agreement|precision|recall|f1|
|---|---|---|---|---|---|---|---|---|---|
|bp-training.csv|47691|47691|0|1|1|1|1|1|1|
|bp-validation.csv|5252|5252|0|1|1|1|1|1|1|
|bp-test.csv|2392|2392|0|1|1|1|1|1|1|
|cc-training.csv|45309|45309|0|1|1|1|1|1|1|
|cc-validation.csv|4985|4985|0|1|1|1|1|1|1|
|cc-test.csv|1265|1265|0|1|1|1|1|1|1|
|mf-training.csv|32421|32421|0|1|1|1|1|1|1|
|mf-validation.csv|3587|3587|0|1|1|1|1|1|1|
|mf-test.csv|1137|1137|0|1|1|1|1|1|1|

## Pickle Summary

|file|status|generated_exists|reference_exists|generated_shape|reference_shape|protein_jaccard|annotation_or_term_jaccard|
|---|---|---|---|---|---|---|---|
|train_data.pkl|compared|True|True|(66841, 3)|(66841, 3)|1|1|
|test_data.pkl|compared|True|True|(3328, 3)|(3328, 3)|1|1|
|terms.pkl|compared|True|True|(5220, 1)|(5220, 1)||1|
|train_data_train.pkl|compared|True|True|(60156, 3)|(60156, 3)|1|1|
|train_data_valid.pkl|compared|True|True|(6685, 4)|(6685, 4)|1|1|

## Known Policy Gaps To Consider

- CAFA3 README includes TAS and IC while some public benchmark code keeps only EXP, IDA, IPI, IMP, IGI, IEP.
- CGD Candida backfill removal may be undocumented outside official materials.
- MF protein-binding-only removal for GO:0005515 may not be implemented by the public Python path.
- GO ontology release choice and obsolete/alternate GO handling may affect propagated labels.
- Reference CSVs may include preprocessing decisions made outside the public benchmark-construction code.
