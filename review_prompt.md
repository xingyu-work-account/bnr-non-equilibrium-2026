# Independent review prompt: non-equilibrium analysis pipeline

## 1. Your role

You are an independent reviewer of a scientific data-analysis pipeline built from `implementation_prompt.md` in this repository. You did not write the code, and you must not trust the implementer's report. Every conclusion you state must rest on evidence you produced in this session.

**Constraints.**

- Do not edit any file in the repository.
- Put any scratch scripts you need in a temporary folder outside the repository, and run them with `uv run --project <repo> python <script>`, so they use the project environment without changing it.
- Never modify `data/` or `outputs/` by hand.

## 2. Scope of review

1. **The complete cumulative change.**
    - If the repository uses git, review `git status` and `git diff` against the state before implementation, including untracked files.
    - Otherwise, list every file outside `data/` and `.venv/` and read all source, test and configuration files in full.
2. **The environment.**
    - `pyproject.toml` and `uv.lock` are consistent.
    - Dependencies were added with `uv`, and no `pip` artifacts or hand-edited dependency lists are present.
    - The build-system change, if any, is correct, and `uv sync` reproduces the environment.
3. **The implementer's report.** Treat every claim in it as unverified until you reproduce it.

## 3. What to verify

1. **Runtime execution.**
    - Run `uv sync`, `uv run ruff check .` and `uv run pytest -q`.
    - Run `rm -rf outputs && uv run python -m bnr_neq run-all` twice, and compare the output tables from the two runs for determinism.
    - Report exact outputs.
2. **Data integrity.**
    - Recompute independently in a scratch script, using only pandas and NumPy on the raw files, without importing `bnr_neq`:
        - row counts;
        - dropped rows;
        - time spans;
        - stamp sizes;
        - the conductivity conversion EC (mS/cm) = −2.89855 × I (A) + 0.520870.
    - Confirm that no row is lost silently, and that reconstructed times are strictly increasing.
    - Confirm that no interpolation crosses a gap longer than twice the native interval.
3. **Experimental facts.** Confirm that the configuration encodes the facts exactly:
    - a 100 mL vessel, open to air and magnetically stirred;
    - 25 mg NH₄Cl at days 8.53, 10.57, 24.59 and 32.39;
    - NH₄⁺ in mg NH₄⁺/L, with a factor of 0.77649 to mg N/L;
    - Pt against a bare Ag reference.

   Confirm that the code reads these values from the configuration and does not hard-code them.
4. **Scientific and numerical correctness.** Recompute independently:
    - the pulse dose: 4.6737 mmol/L, 84.31 mg NH₄⁺/L and 65.46 mg N/L;
    - the Kohlrausch predictions: 0.700 and 0.274 mS/cm;
    - the four NH₄⁺ response fractions;
    - the affinity steps: 4.2, 2.5, 3.1 and 2.4 kJ/mol;
    - the coexistence affinity: about 317 kJ/mol for the documented scenario;
    - pK_a at 25 °C;
    - the air-equilibrium bicarbonate at pH 7.3;
    - the gain-corrected rates after pulses 2 and 4.

   Check the units and signs in every equation in `constants.py` and `thermo.py` against `implementation_prompt.md`, Section 6.
5. **Silent fallbacks.**
    - Search for default arguments, `fillna`, `try`/`except` blocks, `errors="ignore"`, and any code path that substitutes a value for a missing measurement.
    - Every missing input must produce an exception or an explicit `not computed` entry.
    - Scenario values must be labeled as scenarios in every table and figure.
6. **Estimator validity.**
    - Rerun the Markov-chain validation independently with different random seeds and (p, q) pairs, including p = q.
    - Confirm that the irreversibility estimators recover (p − q) ln(p/q) and give zero for reversible surrogates.
    - Confirm that ties in pH are handled as specified.
7. **Test adequacy.**
    - The tests exercise the production functions.
    - They include negative and boundary cases: unit errors, missing inputs, non-monotonic time and impossible chemistry.
    - No test was weakened, skipped or made tautological.
    - Name any untested production path.
8. **Over-claiming.**
    - Confirm that the report labels conclusions as GO, GO WITH DOCUMENTED NON-BLOCKING DEBT or NO-GO, as in the implementation prompt (Section 8, M8).
    - Confirm that no result relies on conductivity as solution conductivity, on absolute ORP, or on energy without stirrer power.
    - Confirm that NH₄⁺ relaxation times are not presented as process times before the electrode step response is known.
9. **API and schema.** Check:
    - the command-line commands;
    - the output file names and column names against the prompt;
    - that units appear in column names or headers.

## 4. Severity and decision

Classify every finding by severity:

- **P0:** invalidates a conclusion or loses data. Examples: a unit error, silent data loss, a wrong equation, a rate computed across a pulse, a second-law violation, or a non-reproducible output.
- **P1:** undermines a specific claim. Examples: an unvalidated estimator, missing uncertainty, a missing negative test, an unlabeled scenario value, or train and test leakage.
- **P2:** documentation, naming or style.

**Decision rule.**

- Any unresolved P0 or P1 finding, failed validation command, unverified acceptance criterion, unresolved scientific or data-integrity risk, or broken reproducibility requires **NO-GO**.
- Otherwise, issue **GO**, or **GO WITH DOCUMENTED NON-BLOCKING DEBT** if P2 items remain.

## 5. Report format

1. **Decision.** GO, GO WITH DOCUMENTED NON-BLOCKING DEBT, or NO-GO, in one sentence.
2. **Commands run.** Each command with its exact output.
3. **Independent recomputation.** A table of each quantity, your value, the pipeline's value and the agreement.
4. **Findings.** Ordered P0, then P1, then P2. Each finding gives:
    - the file and line;
    - the evidence;
    - a concrete failure scenario;
    - the required fix.
5. **Acceptance criteria.** Each milestone gate in the implementation prompt, marked verified or not verified, with evidence.
6. **Remediation prompt.** If the decision is NO-GO, a remediation prompt for the implementer that lists only the findings to fix, the tests to add, and the commands that must pass. It must be written so that a fresh agent can execute it without this review.
