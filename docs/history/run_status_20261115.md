# Run status: November 15, 2026, approximately 6:46 a.m. Pacific

**Update, approximately 10:05 a.m. PST:** public-login SSH access has recovered. Live scheduler queries now confirm the generation and AF3 jobs remain queued, as detailed below. The early-morning access failure is retained as history rather than a description of the current blocker.

## Confirmed progress

The new glycan integration preflight **1693262** completed successfully in 1m44s on GPU node g8. Its synchronized log records `passed: true`, `full_forward_backward: true`, finite gradients, a repulsive-gradient check and a zero distant-control penalty. The actual integration-check timestamp is **November 14 at 9:28 p.m. PST** (`2026-11-15T05:28:12.100968+00:00`). This verifies the integration test, not a successful designed binder.

The last successful scheduler snapshot is **November 15 at 2:34:58 a.m. PST** (`2026-11-15T10:34:58.185781+00:00`). It records input validation and glycan preflight as completed; Ig1/Ig2 generation, their downstream quality/triage jobs and all AF3 jobs were pending. The last HTML build was approximately 2:35 a.m. and still contained the eight historical variants from four Ig1 backbones, zero qualified models and two comparators awaiting predictions.

## Current access failure

Current SSH attempts fail with:

```text
Access denied by pam_slurm_adopt: you have no active jobs on this node
```

The error was reproduced through the normal `hoffman2.oarc.ucla.edu` endpoint and all four DNS-listed login nodes: login1–login4 (`192.154.2.205`–`192.154.2.208`). Alternate-node checks retained strict verification against the existing cluster host key. No credential, known-host entry or SSH configuration was changed.

The local monitor process **69267** remains alive, but its recent polls fail with SSH exit 255. Therefore **the current remote job states are unknown**. The cached pending states must not be presented as a fresh scheduler observation, and the SSH message does not establish whether the batch jobs have run, completed or failed since the last successful check.

Slurm documents `pam_slurm_adopt` as an SSH access-control module for nodes with running allocations ([official documentation](https://slurm.schedmd.com/pam_slurm_adopt.html)). UCLA lists `hoffman2.oarc.ucla.edu` as the public login endpoint ([official connection documentation](https://www.hoffman2.idre.ucla.edu/Using-H2/Connecting/Connecting.html)). Seeing this denial on all public login nodes warrants an administrator check; the cause has not been confirmed. It is not a model-inference traceback.

No job was canceled, resubmitted or reconfigured during this review. Existing local work and frozen campaigns were preserved. The watcher continues its existing retries; a new watcher was not started.

## Support-ready description

User `sujit009` cannot SSH to the production Slurm login endpoint `hoffman2.oarc.ucla.edu` or login1–login4 as of November 15 around 6:46 a.m. PST. Every tested endpoint returns `Access denied by pam_slurm_adopt: you have no active jobs on this node`. These are public login endpoints, not requested direct connections to GPU compute nodes. Please investigate login-node PAM/access policy and confirm the status of the existing jobs below.

Account: `pi_pallard`; partition: `campus`; QoS: `campus24`. Project: `/u/scratch/s/sujit009/c288-mertk-binder/runs/round_20261114_ig_guard_01`.

- Passed glycan preflight: **1693262**.
- Ig1 generation: **1693263**; Ig2 generation: **1693266**.
- Conditional parent AF3: **1693273**.
- Comparator AF3, TULP2: **1693284**; mCherry: **1693312**.
- Other tracked jobs: **1693261**, **1693265**, **1693268**, **1693269**, **1693319**.

Last accessible scheduler snapshot: November 15 at 2:34:58 a.m. PST. The local monitor is still retrying, but no current Slurm status can be retrieved. This description is a draft for [Hoffman2 support](https://support.idre.ucla.edu/); it has not been sent.

## Recovery and live queue check, 10:05 a.m. PST

The standard `hoffman2.oarc.ucla.edu` endpoint accepted SSH again. Live `sacct`, `squeue` and `scontrol` checks confirm:

| Stage | Job | Current state |
|---|---:|---|
| CPU input validation | 1693261 | COMPLETED, exit 0:0 |
| Glycan GPU integration preflight | 1693262 | COMPLETED, exit 0:0, 1m44s |
| Ig1 generation | 1693263 | PENDING, ReqNodeNotAvail |
| Ig2 generation | 1693266 | PENDING, ReqNodeNotAvail |
| AF3 TULP2 comparator | 1693284 | PENDING, ReqNodeNotAvail |
| AF3 mCherry comparator | 1693312 | PENDING, ReqNodeNotAvail |
| Quality, triage, conditional parent AF3 and final report | 1693265, 1693268, 1693269, 1693273, 1693319 | PENDING, Dependency |

The two generation jobs have no remaining scheduler dependency. They became eligible after the successful preflight at approximately 9:29 p.m. on November 14, but have received zero GPU runtime. Both accept A100/H200/V100 and request one GPU, eight CPUs and 64 GB host memory. Comparator jobs request one GPU with 80/141 GB GPU-memory features and 96 GB host memory. The configured account, partition and QoS remain `pi_pallard`, `campus` and `campus24`.

The live reason lists unavailable nodes; Slurm reports `StartTime=Unknown` and `squeue` provides `N/A` for estimated start. `sinfo -R` includes GPU reset/hardware and node-response issues on some nodes, but does not establish that every eligible GPU is down or explain the duration of the wait. No application crash or canceled generation job appears in the current suite.

The existing watcher has resumed successful synchronization. Its 10:04 a.m. snapshot and the rebuilt HTML report the same eight historical variants/four backbones, zero qualified models and two comparators awaiting predictions, with no build errors. There are no new generation/AF3 logs or structures yet. No job or resource request was changed during this check.
