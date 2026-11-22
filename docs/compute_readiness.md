# Compute readiness snapshot

Checked November 10, 2026. This records available hardware and prerequisites, not a completed inference benchmark.

## Local machine

The inspected MacBook Pro has an Apple M5 Pro, 18 CPU cores, a 20-core Apple GPU, 48 GB unified memory and approximately 70 GB free disk space. No NVIDIA/CUDA runtime is present. Unified Apple GPU memory does not provide CUDA compatibility.

| Stage | Local feasibility | Current status |
|---|---|---|
| Reference preparation, structural analysis and plots | Supported by the lightweight project environment | Experimental 37NT source and exported target hashes verified |
| BindCraft generation | Upstream requires NVIDIA/CUDA; Apple GPU is not a supported backend | Pinned source and prepared input exist; no generation launched |
| RFD4-Proteina generation | Qualified runtime is Linux with NVIDIA/CUDA; its CPU image supports preprocessing/metrics | Pinned source on Mac; Hoffman2 runtime and weights installed, generation queued |
| LigandMPNN refinement | Upstream explicitly falls back to CPU | Needs installation and actual candidate complexes before refinement |
| AF3 assessment | Experimental CPU/Apple GPU route exists, distinct from its validated NVIDIA path | Requires authorized parameters, a separate runtime and protocol/backend changes; not benchmarked locally |

For this project, the practical route is to run the generators on Hoffman2 and keep input preparation, result inspection and visualizations on the Mac. Local AF3 feasibility must be measured on a small case; no runtime or numerical-equivalence claim is established. The local disk does not accommodate AF3's full genetic database installation, so a local assessment would require prepared MSAs or an explicitly documented single-sequence protocol.

## Hoffman2 preflight

SSH access to the production Slurm login node was verified, including its ED25519 fingerprint against UCLA's published fingerprint. Account inspection reported access to campus GPUs, with a two-GPU user limit. A100 hardware was listed. No jobs were submitted; no compute allocation was consumed. The old `hoffman2` SSH alias points to the legacy Grid Engine cluster, so new runs must explicitly use the production Slurm hostname.

The user selected Hoffman2 for BindCraft, RFD4-Proteina and AF3, with Mac-side analysis and figures. All environments, model downloads and generated runs live under `/u/scratch/s/sujit009/c288-mertk-binder`. After checking the available Slurm association, the user authorized `pi_pallard`. The initial no-submission statement above describes the preflight only; setup and dependent GPU jobs have now been submitted. See [the pilot run record](hoffman2_pilot.md).

AF3 environment setup has completed. Checking the current upstream README revealed an official direct parameter download, and the standard parameters have now been downloaded and validated in scratch. Shared AF3 assessment is queued after generated candidates are prepared. No completed design results are claimed. Initial figures will show native complex structures, receptor contact positions and receptor deformation relative to 37NT. Shared AF3 confidence and pose-consistency comparisons require the same AF3 protocol for both generator arms.

## Evidence

- [RFD4-Proteina installation and CPU-container scope](https://rfd4-proteina.org/latest/installation.html).
- [BindCraft requirements](https://github.com/martinpacesa/BindCraft#installation).
- [LigandMPNN device selection](https://github.com/dauparas/LigandMPNN/blob/26ec57ac976ade5379920dbd43c7f97a91cf82de/run.py#L38).
- [AF3 CPU and experimental Apple GPU installation](https://github.com/google-deepmind/alphafold3/blob/28c2d290559fe07387a12f471ad49f0d36206327/docs/installation.md#running-alphafold-3-without-a-gpu-cpu-only).
- [Hoffman2 migration notice](https://www.hoffman2.idre.ucla.edu/News-archive/News.html) and [Slurm resource guidance](https://www.hoffman2.idre.ucla.edu/Using-H2/Computing/Computing-Slurm.html).
