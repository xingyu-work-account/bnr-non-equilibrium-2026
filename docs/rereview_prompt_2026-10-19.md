# Independent re-review prompt for the 2026-10-19 remediation

You are an independent reviewer who did not implement the change. Inspect the complete cumulative diff since commit `eb74b94`, not only the latest commit, and reconcile every claim in the implementer's report against the code and the executed outputs.

## Procedure

1. Run `uv run ruff check`, `uv run pytest -q` and `uv run python -m bnr_neq run-all` yourself and record the full results; a failed command is NO-GO.
2. For each finding R1 to R12 in `docs/review_2026-10-19.md` and each work item W1 to W8 in `docs/remediation_prompt_2026-10-19.md`, state whether it is resolved, partly resolved or unresolved, with file and line evidence.
3. Trace one pulse end to end: confirm that its rate, FA, η, oxygen demand and energy metrics derive from the same inverted NH4 concentration, and that the provenance column says so.
4. Recompute by hand the W2 positive-test values (alkalinity bound, air-equilibrium HCO3−, FA) and one aeration-equivalent energy value from its stated SAE scenario, checking units throughout.
5. Check the runtime paths for silent fallbacks: a missing measurement must produce a status string, never a default number, and no scenario value may appear without its label.
6. Check API and schema changes in `outputs/tables`: every renamed or removed column is documented, and downstream readers (report, figures, tests) were updated.
7. Judge test adequacy: positive, negative, boundary, regression, integration and end-to-end tests exist for each new module and exercise the production path rather than a copy of it.
8. Judge scientific correctness: the process-state rules were fixed before evaluating associations, block bootstrap respects autocorrelation, the alkalinity bound is not contradicted by the reported nitrification extent, and the report does not present consistency checks as findings or associations as causes.

## Decision

Classify each remaining issue as P0, P1 or P2. Any unresolved P0 or P1, failed command, unverified acceptance criterion, unresolved scientific or data-integrity risk, or broken compatibility gives NO-GO. Otherwise decide GO, or GO WITH DOCUMENTED NON-BLOCKING DEBT listing each P2 item, and justify the decision with evidence.
