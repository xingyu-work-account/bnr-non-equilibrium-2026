# Independent re-review prompt for remediation pass 6 and W6/W7 (2026-11-05)

You are an independent reviewer who did not implement the change. The previous re-review (`docs/rereview_2026-10-25_pass5.md`) decided GO WITH DOCUMENTED NON-BLOCKING DEBT, and pass 5 was merged into `main` as `0f7abfc`. The implementer then made five commits on branch `remediation/pass6-2026-10-26`, from `2dd5a2b` (whose code and outputs equal `0f7abfc`):

| Commit | Content | Claimed effect on outputs |
| --- | --- | --- |
| `f1744fc` | G1/G2: `docs/states_registration.md`, the SHA-256 registration of the W3 rules, the bootstrap null-coverage check and its decision rule, committed **before** they were run | none (code and tests only) |
| `ba394df` | G3, G4, G6–G8: SCENARIO/ASSUMPTION labels and NH4 observation-model provenance columns; the G2 results | new columns and one new table; pre-existing values unchanged except the `thresholds` text of `states_rules.csv`; fig4 panel-c title only |
| `d6d9e5d` | G5: the aeration-efficiency (AE) source rule, written down before any edit | none (log only) |
| `da49e95` | G5: field AE range 0.7 / 1.6 / 3.5 kg O2/kWh from Crawford and Sandino (2010), WERF Table 2-8 | only the 108 AE rows of `pulses.csv`, `pulses.md`, `analysis_report.md`, fig1 and fig8 |
| `0178cd8` | W6/W7: Figures 1, 4 and 8; 29 single-axes figures plus `no_text/` copies; report reordered; consistency-check labels; rounding | no table changed; report order, rounding and figure files only |

Decide whether these commits (a) close G1 to G8 of the pass-6 audit, (b) complete W6 and W7 of `docs/remediation_prompt_2026-10-19.md`, and (c) change nothing that was accepted in passes 1 to 5.

Treat everything in `docs/progress.md`, sections "Remediation pass 6", "Remediation pass 6, G5" and "W6 and W7", as claims to verify, not as facts. Use `git diff 2dd5a2b 0178cd8` for the whole change, and the per-commit diffs for each guard below.

**There is no separate remediation prompt for pass 6.** Its scope (G1–G8) comes from an audit of W3 to W5 recorded only in `docs/progress.md`. Part of your job is to judge whether that audit was complete.

Reference documents:
- `docs/remediation_prompt_2026-10-19.md`: work items W1 to W8, the invariants, tests, acceptance criteria and failure conditions. W3 to W7 are the ones in scope.
- `docs/rereview_2026-10-25_pass5.md`: the accepted state, open debt P2-1 and P2-2, and the rulings of §6. Do not reopen rulings already made there.
- `docs/states_registration.md`: the registered hash and the stated history of the W3 rules.
- `docs/measurement_request_2026-10-26.md`: context only; not under review.

## Procedure

1. **Run the commands yourself** on a clean checkout of `0178cd8`, and record the full results. A failed command is NO-GO. The implementer's runs were long on a shared machine, so plan for several hours, and do not run pytest and run-all at the same time.
   - `uv sync` and `uv run ruff check`.
   - `rm -rf outputs && uv run python -m bnr_neq run-all`, with a worker count different from the implementer's default (all cores; set `BNR_NEQ_WORKERS`). The implementer reports exit 0 after 3:08:32, with stage lines unchanged from pass 5.
   - `uv run pytest -q -rs`, after run-all. The implementer reports `491 passed`, with none skipped. **That run came before a ruff rename in `tests/test_figures.py` (`out, seen` → `_, seen`), was not repeated, and recorded no code hashes.** Your run on `0178cd8` is the first full test run of the committed code. Record the SHA-256 of every file under `src/` and `tests/`, plus `config/experiment.toml`, `docs/states_registration.md`, `pyproject.toml` and `uv.lock`, before run-all and after pytest.
   - Diff the regenerated `outputs/` against a copy of the committed outputs taken before the run. Everything except `.DS_Store` must be byte-identical.

2. **Guards: nothing accepted may change.** Check each one on your own regenerated outputs, and against the committed outputs of the stated baseline.
   - **G1–G4, G6–G8 against `2dd5a2b`.** In every table they touched, pre-existing columns must be identical as raw CSV text, with row counts and column order unchanged. The only exception is the `thresholds` text of `states_rules.csv` (G3). The only new columns are the ones listed in `progress.md`, plus the new table `states_bootstrap_null_coverage.csv`. `processed/` must be byte-identical. Check `git diff 2dd5a2b ba394df -- src` for any change to a fit, bound, start, selection rule, state rule or primary computation.
   - **G5 against `ba394df`.** The 244 rows of `pulses.csv` without an AE suffix must be byte-identical. For the 108 AE rows, recompute the scaling yourself for at least one row of each of the five metric kinds: AE_old/AE_new for energy and AE_new/AE_old for the Gibbs-to-aeration ratios, to 1e-9 relative. Confirm that no other table changed.
   - **W6/W7 against `da49e95`.** No table may change. `0178cd8` also touches `dynamics.py`, `model.py`, `pulses.py`, `segments.py` and `validation.py`. Confirm, function by function, that those edits only move or delete plotting code (into `figures.py`) or rename docstring cross-references, and that no computation changed. Confirm that the deleted figure stems `fig{1..8}_*` are no longer referenced anywhere in the code, tests or reports.

3. **Close each gap, with file and line evidence.**
   - **G1.** Recompute `states.registration_sha256()` and compare it with `docs/states_registration.md`. Check that the hash input covers what the document says: the rules, thresholds, windows, seeds, resample counts, exclusions, partitions, and the source of every function in `states._REGISTERED_FUNCTIONS`. Make one change of your own to a decision function and confirm the test fails, then revert it. Check with `git log` that `f1744fc` (2026-10-26 16:02 −04:00) contains the null-check code and its 0.90 decision rule, and that no output of the null check appears in any earlier commit. Confirm that every output which reads the rules carries `REVISED AFTER FIRST RESULTS (not pre-registered)`, and that the history in the registration document (the pass-2 `MIN_BLOCKS_PER_GROUP` rule and the pass-3 X2 widening) matches the earlier sections of `progress.md`.
   - **G2.** Reproduce at least one cell of `states_bootstrap_null_coverage.csv` independently, for example pH rising vs not (pooled) at φ = 0.9. Judge whether the null design (the 101-bin least-squares slope of AR(1) noise on the record's own bins and group labels) is a fair null for the real response. Confirm that the report does not read the state-vs-state interval (1.02 to 5.83, excluding 0) as evidence of a difference.
   - **G3, G4.** Check the 14.8 mV value (RT ln10/(4F) at 298.15 K). Count the `assumptions` labels in `pulses.csv` (248 / 28 / 16 / 3 × 12 / 24). For at least one row in each category, check from the code path that the label is correct, including FA from the electrode reading being free of the volume assumption. Confirm that the 24 rows relabelled after the first full run (development note 3) are now correct.
   - **G6, G7.** Confirm that K_NH, K_ALK and the closed alkalinity balance are labelled on every M6 table, and that `model.md` states the scope, "an alkalinity limitation, not an inorganic-carbon limitation".
   - **G8.** Confirm that the seven M6 tables, `nh4_truth_vectors.csv` and `nh4_truth_w2_consistency.csv` carry the observation-model and series columns. Confirm that `check_m6_provenance` and `check_single_observation_model` run before every write, by finding every write site. Delete one provenance column in a scratch copy and confirm that the run raises.
   - **G5.** The source PDF for WERF Table 2-8 is **not in the repository**. `reference/` is untracked and holds only the Gresch et al. (2011) and Desalination and Water Treatment PDFs, plus a secondary summary. Obtain Crawford and Sandino (2010), WERF OWSO4R07e, yourself and check the 10 field rows, the footnote "Includes blower power requirements", and the page number. Recompute the rule: the minimum and maximum of the three diffused-air rows marked with that footnote, and their geometric mean rounded to 2 significant figures. Check that the DTW PDF in `reference/` (file name `S1944398624019593`) is the paper cited (Azis@Aziz et al. 2019, doi:10.5004/dwt.2019.24134). Confirm that the rule in `d6d9e5d` predates the code in `da49e95`, and that the generated DO wording differs from the registered wording only as `progress.md` discloses.
   - **W6.** Check figures 1, 4a–c and 8a–d against the W6 row: the linkage chain, the N_nit(Alk0) bound, rate against FA from W1, the state timeline, η, oxygen utilization and the idle fraction under the AE and kLa SCENARIOs. Check that every number printed in Fig. 1 equals its source table or config value. Open the PNGs and judge whether they are readable at the fixed 25 pt font. Check that "n.c." marks and the conditional marker appear only where `pulses.csv` says so.
   - **W7.** Read the regenerated `analysis_report.md` from start to end. The order must be narrative, then energy and efficiency, then sensor validity, with the stochastic-thermodynamic results in a supplement. **Settle one possible conflict:** W7 asks for sensor validity after energy, but `progress.md` says "Sensor validity" moved from §5 to §3, before process state (§4) and the model (§5). Decide whether the final order meets W7. Confirm that the coexistence affinity, ΔA and σ ≥ 0 are labelled CONSISTENCY CHECK everywhere they appear. Spot-check the rounding: no non-constant value with more than 4 significant figures.

4. **Judge test adequacy.**
   - Check `tests/test_pass6_fixes.py` (31 tests, after the G5 test) and `tests/test_figures.py` (2 tests) against the positive, negative, boundary, regression, integration and end-to-end rows of `docs/remediation_prompt_2026-10-19.md`.
   - For every changed existing test, compare it with `git show 2dd5a2b:<file>`: `test_pass4_fixes.py`, `test_pass5_fixes.py`, `test_report.py`, `test_segments.py` and `test_validation.py`. Decide whether each change only follows appended columns, new section numbers or new figure stems, or whether it weakens an assertion. In particular, the extended column-equality tests must still fail on any change to a pre-existing column.
   - Decide whether `test_figures.py` checks the plotting rules on the saved files, not only on the figure objects.

## Rule explicitly on the implementer's judgment calls

1. **G1, freezing after the fact.** The rules are hashed and labelled "revised after first results" instead of being presented as pre-registered. Is that sufficient, and is the label used consistently?
2. **G2 design and threshold.** The parameters are φ ∈ {0.5, 0.9, 0.95, 0.98}, 200 replicates, 999 resamples, and "adequate" at a coverage of 0.90 or more at every φ. The block rule stays as registered whatever the check finds. Two stage previews ran after `f1744fc` but before the validation run, and only report rendering changed between them.
3. **G4 correction after the first full run.** 24 rows were wrongly labelled volume-free and were fixed together with the unit test that encoded the same error. Is the strengthened test (every row computed from C_inv or the model state names its volume frame) adequate?
4. **G6/G7, no rename and no CO2 term.** The candidate names stay as they are, and the scope is stated in labels and in `model.md`.
5. **G5 source rule.** Only blower-inclusive diffused-air rows are used, with a geometric-mean centre. The shorter DO wording was used. The Stenstrom and Rosso (2010) "not for design" caveat was dropped from the report. Baeten et al. (2021) is left unused because it could not be verified.
6. **W6 single-axes split.** Multi-panel figures became 29 single-axes figures with a `no_text/` copy each. Arial has no weight 600, so it is drawn at 700.
7. **Validation of the final commit.** pytest was not rerun after the ruff rename, and no code hashes were recorded for the W6/W7 run. Your step 1 replaces that evidence. Rule on whether the omission matters beyond this review.
8. **Status of P2-1 and P2-2 from pass 5.** Say whether either was touched or is still open. Pass 6 did not claim to address them.

## Out of scope, but report if found

- `docs/reports/BNR_non_equilibrium_full_report_2026-10-28.{md,html}` is untracked. It claims to add no new computation. If time allows, spot-check that five of its headline numbers equal `analysis_report.md` and the tables, and note any that do not. Its status does not affect the decision.
- `src/bnr_non_equilibrium_2026/` and `reference/` are untracked. Say whether either should be committed, for example the G5 source material, for the record to be reproducible.

## Decision

Classify each remaining issue as P0, P1 or P2. Any of the following gives NO-GO:
- an unresolved P0 or P1;
- a failed or non-reproducible command, or a test that fails at `0178cd8`;
- a pre-existing value, figure body or `processed/` file that changed beyond what is disclosed and explained above;
- a registration hash that does not match the code, or a null-check output older than `f1744fc`;
- a G5 range that does not follow from WERF Table 2-8 under the registered rule;
- an unverified acceptance criterion, or a met failure condition, of `docs/remediation_prompt_2026-10-19.md` for W3 to W7, in particular a SCENARIO value presented as measured or an output that mixes observation models;
- a rejected judgment call rated P0 or P1;
- a weakened test.

Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT, listing each P2 item. If the decision is GO or GO WITH DEBT, say whether `remediation/pass6-2026-10-26` is ready to merge into `main`, and what must be committed or pushed first. Write the review to `docs/rereview_<date of your run>_pass6.md`, using the structure of `docs/rereview_2026-10-25_pass5.md`: Decision, Commands, Guards, gap-by-gap closure (G1–G8, W6, W7), Scientific notes, Tests, Rulings, Remaining issues by severity, Readiness to merge, and What is sound and should be kept.
