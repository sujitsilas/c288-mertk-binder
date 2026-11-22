# Stalk alignment sensitivity

## Audit the geometry screen

The screen slides and rotates each predicted complex to match the 14 Cα atoms at MerTK residues 485–498. Those points extend mainly along one axis. A small change in their coordinates can rotate an off-axis binder appreciably, even when the fitted Cα atoms remain close.

The diagnostic applies the same 200 draws of independent, zero-mean Gaussian noise to each conformer’s Cα coordinates. The standard deviation is 0.5 Å per Cartesian coordinate. Each noisy cloud is aligned to its original cloud with the production rigid fit. Two probes lie 20 Å from the centroid along the minor principal axes. The metric is the larger displacement of the two probes per draw. The random seed is 2436. These settings define a numerical stress test. The displacements are neither experimentally measured errors nor confidence intervals for individual designs.

| Stalk conformer | Singular values of centered coordinates (Å) | 95th-percentile probe displacement (Å) | Largest observed displacement (Å) |
| --- | --- | ---: | ---: |
| Base | 48.735, 6.696, 3.272 | 2.623 | 4.093 |
| P490 ψ +25° | 49.220, 5.113, 3.223 | 3.406 | 4.714 |
| P490 ψ −25° | 47.200, 9.414, 3.389 | 2.112 | 2.647 |

The probe sensitivity is comparable to the 2.0 Å protein and 2.5 Å glycan distance thresholds. Interpret the recorded 0 of 50 cross-conformer passes as the result of these geometric checks on three chosen models. The test does not establish whether any candidate binds a native receptor. Each future geometry result includes the coordinate spread of both aligned shapes.

The median fraction of raw predicted contacts on residues 485–498 was 0.264. Forty-six of 50 poses had fewer than half of their contacts on the stalk, with 48 failing the combined contact test for insufficient cleavage-site engagement or predominantly off-stalk contacts. Of 150 pose–conformer combinations, 33 fell below the assumed 5 Å membrane clearance margin. Nine crossed the estimated membrane plane. Both counts depend on the modeled membrane orientation.

## Reproduce the audit

```bash
.venv/bin/python scripts/audit_stalk_fit.py \
  outputs/hpc/round_20261117_protector_01/references \
  outputs/hpc/round_20261120_mertk_pose_followup_01/stalk_geometry.json
```

[Audit settings, results, and source hashes](../assets/figures/stalk_alignment_sensitivity.json)

## Contact reporting

`target_contacts` counts native-numbered residues within 4.5 Å of the partner in the original predicted coordinates across the full target construct. `reference_contacts` limits those contacts to residues represented in the overlay reference. Earlier viewers used `contacts` for that reference-limited list. Experimental structure 37NT ends at Gly484, so its reference-limited contacts omit the modeled stalk. Stalk counts use `target_contacts`.
