# Blackbox Task 50 - Förster Transfer Efficiency in the Doubly Static Regime: Frozen Isotropic κ² and a Radial-Gaussian Distance Ensemble

**Domain:** Chemistry - physical chemistry / photochemistry (resonance energy transfer, orientation statistics, conformational ensembles)
**Tier:** Med (floor: 3h, 2 genuine twists, 3 substantive steps)
**Honest counts:** 3 substantive steps (plus a Förster kernel and a trivial ×100/rounding step, not counted), 2 genuine twists

**Repo name:** `chemistry-fret-static-kappa-radial-gaussian-efficiency`
**GitHub description:** Reverse-engineer the ensemble-averaged Förster transfer efficiency of a donor-acceptor pair whose dipole orientations and separations are both frozen on the transfer time scale. Requires the static isotropic κ² distribution, an r²-weighted distance ensemble, and a first-principles Förster radius. The textbook ⟨κ²⟩ = 2/3 answer misses by up to 19 percentage points.

---

## Task description (for reviewers)

The hidden function maps four positive reals to a percentage. It is the efficiency of a distance-dependent photophysical process, averaged over a population of molecules. The oracle builds the characteristic length from first-principles constants. It then averages the per-molecule efficiency, not a rate constant, over two frozen sources of heterogeneity: one in molecular orientation and one in separation. A solver must identify the process, find the non-textbook short-range law (100 − output ∝ r³ rather than r⁶), deduce which orientation statistics produce it, and find the Jacobian-weighted form of the distance ensemble. The target is 0.001 percentage points. Expected solver knowledge: graduate photophysics (Förster theory, orientation-factor statistics) and polymer conformational statistics.

---

## 1. Domain & Algorithm

Physical chemistry / photochemistry. The quantity is the **population-averaged Förster resonance energy transfer (FRET) efficiency, in percent**, for a donor-acceptor pair. Both its transition-dipole orientations and its D-A distance are *static*: they do not change during the donor excited-state lifetime. Each molecule transfers with its own efficiency E(κ², r) = κ²X/(κ²X + r⁶). The observable, e.g. steady-state donor quenching 1 − F_DA/F_D, is the average of E over the joint ensemble.

## 2. Core Method

- **Förster radius from first principles** (Förster 1948). R₀⁶ = 9000 ln10 κ² Q_D J / (128 π⁵ N_A n⁴), with J in M⁻¹ cm⁻¹ nm⁴ converted to cm⁴ (10⁻²⁸) and the result to Å⁶ (10⁴⁸). This gives R₀⁶ = κ² X with X = 8.785×10⁻⁵ n⁻⁴ Q_D J Å⁶ (Lakowicz eq 13.5: 8.79×10⁻⁵). Refractive index n = 1.4, the conventional value for biomolecules in aqueous solution.
- **Per-molecule efficiency.** E = k_T/(k_T + 1/τ_D) = 1/(1 + r⁶/(κ²X)).
- **Static isotropic orientation factor** (Dale, Eisinger & Blumberg 1979; van der Meer 2002; Loura 2012 eq 1). With D and A dipoles independently isotropic and frozen, κ² ∈ [0, 4] has the density
  p(κ²) = ln(2+√3) / (2√(3κ²)) for 0 ≤ κ² ≤ 1, and p(κ²) = ln[(2+√3)/(√κ² + √(κ²−1))] / (2√(3κ²)) for 1 ≤ κ² ≤ 4, with ⟨κ²⟩ = 2/3.
- **Radial-Gaussian distance ensemble** (Haas, Wilchek, Katchalski-Katzir & Steinberg 1975; in modern form Koren et al. 2023 eq 4). P(r) = K r² exp[−b(r − a)²], normalised so that ∫₀^∞ P dr = 1. The oracle uses b = 1/(2σ²).

## 3. Twists (2 genuine)

**T1 - Static κ² averaging of the efficiency (not of κ²).** The textbook practice (Lakowicz; nearly all R₀ tables) inserts ⟨κ²⟩ = 2/3 into R₀. That is exact only in the dynamic isotropic limit, where orientations randomise faster than transfer. In the static regime each molecule keeps its own κ², and the observable averages the *non-linear* E(κ²). Because p(κ²) diverges as (κ²)^−½ at 0, a finite fraction of molecules barely transfer however close they are. The result:
- 100 − ⟨E⟩ ∝ (r/R₀)³ at short range, rather than (r/R₀)⁶.
- At r = R₀(2/3) the efficiency is 37.9 %, not 50 %.
- Up to 19 percentage points below the textbook value across the domain.

In the far regime E is linear in κ², so the static and dynamic averages coincide to leading order. T1 fades there (control cases).

Implementation: substituting u = |κ| and integrating by parts across the logarithmic kink at κ² = 1 gives a smooth one-dimensional integral,
⟨E⟩_κ(r) = (1/√3) ∫₀^{arccosh 2} [cosh t − √c · arctan(cosh t/√c)] dt, with c = r⁶/X,
evaluated by 48-point Gauss-Legendre. It is exact to machine precision and was checked against a 200 000-sample Monte Carlo over random dipole pairs.

**T2 - The distance ensemble is a *radial* Gaussian with an r² Jacobian, clipped at r = 0.** value_a is the centre a of exp[−(r−a)²/2σ²], not the mean distance. The r² factor shifts weight outward (mean ≈ a + 2σ²/a) and empties short distances. A plain Gaussian in r, the other common choice (e.g. Lakowicz ch. 14), puts the effect of σ in the *opposite* direction near r ≈ R₀:
- At (40, 8, 0.5, 10¹⁵), the radial form gives 42.43, the plain Gaussian 49.33, and a single distance 48.60.
- At small a and large σ the difference reaches 28 points.

The average is taken over the κ²-averaged kernel, so T1 and T2 interact. The σ effect is −6.2 points with static κ² and −10.2 points with dynamic κ² at the same point.

Both twists are I/O-isolable:
- value_b = 0 makes T2 bit-exactly neutral (verified against the frozen variant).
- Large value_a/R₀ (small value_c·value_d) makes T1 neutral within tolerance (edge_control_3: 0.0029 both ways).
- Probing value_a at value_b = 0 exposes T1's r³ tail.
- Probing value_b at fixed value_a exposes T2.

## 4. Multi-Step Pipeline

| Step | Computation | Inputs involved |
|---|---|---|
| S1 | Förster prefactor X = 9000 ln10 Q_D J/(128π⁵N_A n⁴) with unit conversion (R₀⁶ = κ²X) | value_c, value_d |
| S2 | Static isotropic κ² average of E(κ², r) (T1): closed-form reduction plus Gauss-Legendre | S1, r |
| S3 | Radial-Gaussian average of S2 over r ∈ [max(0, a−10σ), a+10σ], 16 panels × 32-point GL, normalised by the same rule (T2) | value_a, value_b |
| (S4) | ×100, round to 4 decimals | trivial, not counted |

Each step is separable and substantive. S1 fixes the absolute length scale (the solver must identify the (Q·J)^{1/6} scaling and the constant). S2 and S3 are non-linear averages that cannot be folded into an effective R₀: an effective R₀ fitted at one distance fails at another.

## 5. Input Schema

```json
{
  "value_a": "float [15, 80]",        // centre a of the radial Gaussian, Angstrom
  "value_b": "float [0, 12]",         // width sigma, Angstrom (0 = single distance)
  "value_c": "float [0.05, 1]",       // donor quantum yield Q_D
  "value_d": "float [1e13, 1e16]"     // spectral overlap J, M^-1 cm^-1 nm^4
}
```
R₀(κ² = 2/3) spans 14.0 to 73.1 Å over the domain.

## 6. Output Schema

```json
{ "output": "float, 4 decimals" }     // ensemble FRET efficiency in percent, 0 .. 100
```

## 7. Hardness & Twists

These baselines are measured on the 51 golden cases (`.scratch/baselines.py`, tolerance 0.001):

| Model | Pass | Median abs err | Max abs err |
|---|---|---|---|
| Textbook: single r, κ² = 2/3 | 3/51 | 5.40 | 53.42 |
| Dynamic κ² = 2/3 with full radial distribution (T1 off) | 4/51 | 4.75 | 18.27 |
| Static κ², single r = value_a (T2 off) | 8/51 | 0.25 | 41.46 |
| Static κ², plain Gaussian without r² | 7/51 | 1.05 | 27.63 |
| Dynamic κ², plain Gaussian | 3/51 | 7.83 | 36.73 |
| Full model (example solver, 3 probes after hypothesis) | 51/51 | 8e-6 | - |

Why this is hard:
- **κ² = 2/3 is a reflex.** Almost every worked FRET calculation uses it. The static regime needs the full p(κ²), including its log-kinked second branch on [1, 4], and requires averaging E rather than κ². Using the "static effective κ²" of 0.476 (= ⟨|κ|⟩², verified numerically) quoted in some reviews also fails: that value belongs to a 3D acceptor ensemble, and no single κ² reproduces the r-dependence.
- **The diagnostic is in a tail.** The r³ law shows up only as 100 − output at short distance, where outputs are 95-99.9. A solver must look at log(100 − y) vs log r to see it.
- **The r² Jacobian is invisible without deliberate probing.** Its effect reverses sign relative to the plain-Gaussian intuition near R₀, and it depends on value_a/value_b, not on value_b alone.
- **No linear or additive surrogate works.** The output is sigmoidal in log r, with two different tail exponents, and is smoothed by a distribution whose weight depends on a.
- **Cold read:** four neutral positive reals and an output in 0-100 suggest "a percentage". The (value_c·value_d)^{1/6} length scaling is discoverable and points to Förster transfer. Nothing in the schema reveals the orientation statistics or the form of the distance ensemble.

## 8. Constants & Citations

Every numeric literal in `oracle/implement.py` is whitelisted with its justification in `scripts/verify_task.py` (AST scan).

| Constant | Value | Status | Source |
|---|---|---|---|
| N_A | 6.02214076×10²³ mol⁻¹ | EXACT | SI 2019 defining constant |
| 9000, 128, π⁵, ln 10, n⁴, r⁶ | - | CITED | Förster (1948) R₀ expression, as quoted in Eis & Lakowicz 1993 eq (2) and Wojcik et al. 2018 eq (2) |
| n | 1.4 | CITED | Refractive index for biomolecules in aqueous solution (Wojcik et al. 2018; Lakowicz 2006 §13.2) |
| 10⁻²⁸, 10⁴⁸ | - | EXACT | nm⁴ → cm⁴, cm⁶ → Å⁶ |
| ln(2+√3) = arccosh 2, √3 | - | CITED / EXACT | Isotropic κ² density (Loura 2012 eq 1; van der Meer 2002); |κ|max = 2 |
| r² factor, ∫₀^∞ P = 1 | - | CITED | Radial Gaussian (Koren et al. 2023 eq 4; Haas et al. 1975) |
| 48, 32, 16, ±10σ, Newton 100 / 1e-16, 0.25, 0.5 | - | STRUCTURAL | Fixed quadrature orders and Legendre-root initial guess cos(π(i−¼)/(n+½)) (Abramowitz & Stegun 22.16.6). Converged: max deviation from a 12σ, 8000-interval Simpson reference is below the 5×10⁻⁵ rounding step. |
| 100 | - | EXACT | percent |

### Citations fetched & verified (quote < 15 words, location)

| Item | Verbatim quote | Where fetched |
|---|---|---|
| R₀ prefactor | "R₀⁶ = 9000(ln10)κ²φD0/(128π⁵Nn⁴) ∫…" | Eis & Lakowicz 1993, *Biochemistry* 32, 7981, eq (2), PMC6897574 |
| n = 1.4 | "n is the refractive index assumed to be equal to 1.4" | Wojcik, Solarczyk & Kulakowski, arXiv:1802.04781, §II below eq (2) |
| κ² density | "the analytical result for isotropic dipoles, which is (e.g., [5])" | Loura 2012, *Int. J. Mol. Sci.* 13, 15252, eq (1), PMC3509639* |
| ⟨κ²⟩ | "the average orientation factor ⟨K²⟩ is 2/3" (dynamic isotropic) | Dale, Eisinger & Blumberg 1979, abstract, PMC1328514 |
| Radial Gaussian | "P(r) = Kr² exp −b(r − a)²" | Koren et al., arXiv:2212.01894, eq (4) |
| Normalisation | "K is a normalization factor ensuring that ∫₀^∞ P(r)dr = 1" | same, text below eq (4) |

\*The PMC HTML renders the radicals of eq (1) ambiguously. The implemented form, ln(2+√3)/(2√(3κ²)) and ln[(2+√3)/(√κ²+√(κ²−1))]/(2√(3κ²)), was confirmed three ways: it integrates to 1, it gives ⟨κ²⟩ = 2/3 exactly (analytic check in `.scratch`), and it agrees with a 200 000-sample Monte Carlo over random dipole pairs to within 1×10⁻³ (`verify_task.py`).

**Citation-code match:**
- `FOERSTER_NUM, FOERSTER_DEN = 9000.0, 128.0`, `N_A = 6.02214076e23`, `N_REFR = 1.4`.
- `T_MAX = acosh(2)` (= ln(2+√3)); `/ SQRT3` for the density.
- `r * r * exp(-((r-a)/sigma)**2/2)` for the radial Gaussian.

`verify_task.py` asserts these anchors:
- prefactor·n⁴ = 8.785×10⁻⁵ (Lakowicz 8.79×10⁻⁵);
- R₀ within 0.2 % of Lakowicz's 0.211(κ²n⁻⁴QJ)^{1/6};
- far limit gives ⟨κ²⟩ = 2/3;
- short-range tail exponent is 3.000.

**Honest notes (CONCERN-level, documented):**
- **Doubly static is a modelling regime.** Real dyes are often partly mobile (anisotropy between 0 and 0.4). The oracle is the rigorous static isotropic limit, not a claim about a specific dye pair. That regime is realistic for rigid labels, frozen glasses, and membrane-embedded probes.
- **Independence of κ² and r.** The joint ensemble is taken as a product measure, the standard assumption when no structural model couples them.
- **Unphysical short r.** The radial Gaussian is clipped at r = 0, as in the cited normalisation, even though real D-A contact distances are ≳ 5 Å. This matters only in the value_a = 15, value_b = 12 corner (edge_boundary_2).
- **Secondary sources.** The Förster (1948) and Haas (1975) originals are scanned or paywalled. Their formulas were verified in open papers that quote them verbatim.

## 9. Edge Cases

| id | value_a | value_b | value_c | value_d | output | Rationale |
|---|---|---|---|---|---|---|
| edge_control_1 | 60 | 0 | 0.5 | 1e14 | 1.5788 | r = 2R₀: near textbook (1.6076) |
| edge_control_2 | 80 | 0 | 0.3 | 1e14 | 0.1738 | r = 2.9R₀: T1 neutral to 4e-4 |
| edge_control_3 | 80 | 0 | 0.05 | 1e13 | 0.0029 | r = 5.7R₀: both twists neutral |
| edge_discrim_1 | 40 | 0 | 0.5 | 1e15 | 48.6019 | T1: textbook gives 65.05 |
| edge_discrim_2 | 30 | 0 | 0.5 | 1e15 | 73.6412 | T1: textbook gives 91.27 |
| edge_discrim_3 | 40 | 8 | 0.5 | 1e15 | 42.4299 | T2: radial lowers E; plain Gaussian gives 49.33 |
| edge_discrim_4 | 80 | 12 | 1 | 1e16 | 27.0273 | T1 × T2 (dynamic + radial: 34.30) |
| edge_discrim_5 | 15 | 12 | 0.3 | 1e14 | 44.1509 | r² Jacobian: plain 71.62, single r 85.61 |
| edge_boundary_1 | 15 | 0 | 1 | 1e16 | 99.1601 | r = 0.2R₀: r³ tail (textbook 99.99) |
| edge_boundary_2 | 15 | 12 | 0.05 | 1e13 | 9.3930 | Smallest R₀, widest ensemble, clipped at r = 0 |
| edge_boundary_3 | 40 | 0.5 | 0.5 | 1e15 | 48.5738 | σ → 0 continuity with 48.6019 |

The golden file adds 40 seeded random cases (seed 50), for 51 cases in total.

---

## GATE 1: Screening

| Test | Result | Notes |
|---|---|---|
| Algebraic collapse | PASS | The output is a double non-linear average. It is not additive or separable in any transform of the inputs. value_c and value_d enter only as a product (a genuine physical degeneracy, kept deliberately so that one probe dimension carries the length scale). |
| Domain recall | PASS (CONCERN) | An expert will likely say "FRET efficiency". The static κ² average of E, the radial r² ensemble, and their composition are not implied by the schema, and no single published method is this pipeline. |
| Twist survives | PASS | T1 changes an asymptotic exponent (6 → 3), which no reparametrisation removes. T2 reverses the sign of the σ effect relative to the plain Gaussian. |
| Genuine twist | PASS (2) | T1 contradicts the universal 2/3 practice. T2 departs from the plain-Gaussian textbook distribution and from the single-distance assumption. The first-principles R₀ prefactor is rigour, not a twist. |
| Two-trap | PASS | value_b is independently probeable (value_b = 0 switches T2 off exactly). T1 is exposed through value_a at value_b = 0. |
| Tier fit | PASS | 3 steps ≥ 3; 2 twists ≥ 2 (Med). |
| PhD authenticity | PASS | Distance distributions from static-regime FRET, and the κ² problem, are active research issues in smFRET and TR-FRET structural biology. |

**Verdict: PASS** (one documented CONCERN on domain recall of the FRET family)

## GATE 2: Verification

| Check | Result |
|---|---|
| Twist integrity (freeze-one tests) | PASS. Over 120 random points, T1 max effect 19.1 points; T2 max effect 39.2 (vs single r) and 28.0 (vs plain Gaussian). The twists interact: −6.17 vs −10.16. |
| I/O isolation | PASS. value_b = 0 ⇒ T2 bit-exact neutral. Far regime ⇒ T1 within tolerance. |
| No invented constant | PASS. The AST scan finds only whitelisted literals. |
| Recall-reconstructability | PASS. The example solver recovers the function exactly (51/51) with 3 probes once the family is hypothesised. Hypothesis-forming probing (tail exponents, product degeneracy, σ response) fits well within the 200-probe budget. |
| Citations fetched & verified | PASS. Six quotes fetched (table above). The rendering ambiguity in eq (1) was resolved by an analytic check and Monte Carlo. |
| Citation-code match | PASS |
| Output/schema neutralisation | PASS: key "output", float, value_a to value_d |
| Tier floor | PASS: 3 steps / 2 twists vs Med 3 / 2 |
| Edge cases ≥ 5 | PASS: 11 total (3 control, 5 discriminating, 3 boundary) |
| PhD-level QA | PASS. Requires orientation-factor statistics, conformational-ensemble modelling, and Förster theory together. |

**Verdict: PASS**

---

## Files

```
task-50/
├── oracle/implement.py          reference oracle (hidden)
├── solution/example_solver.py   probe + fit solver (independent Simpson numerics, 51/51)
├── grader/grade_submission.py   tolerance 0.001 absolute, probe budget 200
├── golden/test_data.json        51 cases (11 edge + 40 random)
├── scripts/helpers.py           shared utilities, twist-freezing variant()
├── scripts/build_golden.py      regenerates golden data
├── scripts/verify_task.py       automated GATE 2 checks (incl. Monte Carlo κ² check)
├── scripts/run_tests.py         verify + oracle self-grade + solver grade
└── PROPOSAL.md
```

## References

- Förster, Th. (1948). Zwischenmolekulare Energiewanderung und Fluoreszenz. *Ann. Phys.* 437, 55-75. doi:10.1002/andp.19484370105
- Dale, R. E., Eisinger, J. & Blumberg, W. E. (1979). The orientational freedom of molecular probes. The orientation factor in intramolecular energy transfer. *Biophys. J.* 26, 161-193. PMC1328514
- van der Meer, B. W. (2002). Kappa-squared: from nuisance to new sense. *Rev. Mol. Biotechnol.* 82, 181-196.
- Loura, L. M. S. (2012). Simple estimation of FRET orientation factor distribution in membranes. *Int. J. Mol. Sci.* 13, 15252-15270. PMC3509639
- Haas, E., Wilchek, M., Katchalski-Katzir, E. & Steinberg, I. Z. (1975). Distribution of end-to-end distances of oligopeptides in solution as estimated by energy transfer. *PNAS* 72, 1807-1811. PMC432635
- Koren, G. et al. (2023). Intramolecular structural heterogeneity altered by long-range contacts in an intrinsically disordered protein. arXiv:2212.01894, eq (4).
- Wojcik, K., Solarczyk, K. & Kulakowski, P. (2018). Measurements on MIMO-FRET nano-networks based on Alexa Fluor dyes. arXiv:1802.04781.
- Eis, P. S. & Lakowicz, J. R. (1993). Time-resolved energy transfer measurements of donor-acceptor distance distributions and intramolecular flexibility of a CCHH zinc finger peptide. *Biochemistry* 32, 7981-7993. PMC6897574 (R₀ eq 2).
- Lakowicz, J. R. (2006). *Principles of Fluorescence Spectroscopy*, 3rd ed. Springer, ch. 13-14 (eq 13.5).
- Abramowitz, M. & Stegun, I. A. (1964). *Handbook of Mathematical Functions*, §22.16 / §25.4.29 (Gauss-Legendre).
