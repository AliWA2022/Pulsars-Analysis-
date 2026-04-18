# Pulsar Analysis Project

A third-year astrophysics lab project covering pulsar data analysis across three objectives: dedispersion, searching for pulsars via FFT, and Crab pulsar timing.

---

## Objective 1 — Dedispersion & DM Optimisation

### Data Structure

The raw data is a **3D array** representing radio intensity as a function of:
- **Frequency channel** (narrow sub-band)
- **Phase bin** (one-phase bin within a folded period)
- **Sub-integration** (time chunk)

The pulsar's approximate period is known, so the data is **folded** around it. The period is divided into **1024 phase bins**, acting as the resolution of the pulse profile. Sub-integrations break the total observation time into ~100 chunks.

| Averaging Type | Description |
|---|---|
| **Fully averaged** | All sub-integrations and frequency channels collapsed into one profile |
| **Time averaged** | Sub-integrations collapsed, frequency channels kept — shows intensity vs frequency vs phase |
| **Frequency averaged** | Frequency channels collapsed, sub-integrations kept |

---

### Dispersion Measure (DM)

Radio pulses arrive at different times across frequency channels due to the **interstellar medium**. The delay for each channel relative to the reference (highest) frequency is:

```
delta_t = K * DM * (1/f_low^2 - 1/f_ref^2)
```

where `f_ref` is the **highest frequency channel** (least affected by dispersion).

---

### Dedispersion Pipeline

1. **`dm_to_bin_shift`** — Converts the DM-induced time delay into a number of phase bins to shift for each frequency channel, using the reference (highest) frequency.
2. **`shift_rows`** — Aligns each frequency channel by rolling it by its corresponding bin shift. Channel 0 shifts by 0, channel 1 by x, channel 2 by 2x, etc.
3. **`dedisperse`** — Combines the above: applies shifts to the time-averaged data, then collapses frequency channels to produce a corrected pulse profile. Only produces a clean profile if the trial DM is correct.

---

### Finding the Optimal DM

The correct DM is the one that produces the **highest Signal-to-Noise Ratio (SNR)** after dedispersion.

#### SNR Scoring Functions

- **`compute_score_at_DM`** — Finds the on-pulse window, creates a Boolean mask over the pulse region, computes baseline and sigma from off-pulse bins, and returns both **peak SNR** and **integrated SNR**.
- **`fixed_window_SNR`** — Uses the peak of the uncorrected pulse profile as a reference to define the window. Can be unreliable if the reference peak is wrong.
- **`floating_window_SNR`** — Uses the integrated SNR from `compute_score_at_DM`. More robust than the fixed window.

#### Coarse DM Scan

The `coarse_dm_scan` function sweeps a DM array, evaluates all SNR functions at each trial DM, plots the results, and identifies the DM at peak SNR.

#### Gaussian Fitting

A Gaussian is fitted to the SNR-vs-DM curve to refine the DM estimate and obtain an uncertainty, using `scipy.optimize.curve_fit` with initial parameter guesses derived from the coarse scan.

---

### Jackknife Uncertainty Estimation

Each sub-integration yields an independent best-fit DM. The **Jackknife** method estimates the DM uncertainty robustly:

1. Iteratively remove one sub-integration at a time (treating it as potential noise).
2. For each iteration, fit a **quadratic curve** around the three points near the SNR peak (peak +/- 1 bin) to find the best DM.
3. Collect the 100 best-DM values into a Jackknife array.
4. Compute the mean and uncertainty from the spread of this array.

This gives a statistically robust uncertainty by accounting for the contribution of each sub-integration independently.

---

## Objective 2 — Searching for Pulsars (FFT-based Period Search)

### Data Structure

The data is a **1D time series** of radio intensity. The time axis is constructed from the sampling time `T_samp`:

```
t_i = i * T_samp,   i = 0, 1, 2, ..., N-1
```

---

### Fourier Transform Period Search

The **RFFT** (`numpy.fft.rfft`) is applied to the intensity time series. For each trial frequency, it asks:

> *"How strongly does the data contain a repeating oscillation at this frequency?"*

- If the trial frequency **matches** the pulsar period, the contributions add **constructively** → large FFT magnitude.
- If the trial frequency is **wrong**, contributions cancel → small magnitude.

The frequency array is computed with `numpy.fft.rfftfreq(N, d=T_samp)`:

```
f_k = k / (N * T_samp),   k = 0, 1, ..., N/2
```

#### Assumptions
- Uniform, known sampling interval `T_samp`
- Observation long enough to resolve the pulsar period
- Stable, regularly repeating pulse signal

> **Note:** The FFT may not give the exact period if it falls **between frequency bins**.

---

### Harmonic Summing

Pulsar pulses are **non-sinusoidal** but periodic, so they decompose into a fundamental frequency `f0` and harmonics `n*f0`. By identifying harmonic peaks and fitting a best-fit line to `f = n*f0`, the **gradient** gives the fundamental frequency `f0`, and hence the period `P = 1/f0`.

The intensity time series is then **folded** at this period to build the final pulse profile.

---

## Objective 3 — Crab Pulsar Timing

### Times of Arrival (TOAs)

A **Time of Arrival (TOA)** is the timestamp of a detected pulse, measured by cross-correlating the observed pulse profile with a **standard template**. TOAs are computed assuming the telescope is at Earth's surface, then corrected to the **Solar System Barycentre (SSB)** to remove Earth's orbital and rotational motion.

---

### Timing Residuals

After barycentre correction, the **timing residual** `t_res` is the difference between the observed TOA and the predicted TOA from a timing model. Plotting `t_res` vs pulse number `N` reveals systematic trends.

#### Unwrapping Residuals

Residuals may contain **phase wraps** (jumps of +/- 1 period). These are corrected by:
- If the residual jumps **up** by ~1 period, subtract one period.
- If it jumps **down** by ~1 period, add one period.

#### Quadratic Fitting & Outlier Rejection

Two rounds of **quadratic fitting** are applied to `t_res` vs `N`:
1. First fit to the full dataset.
2. Remove outliers more than **5 sigma** from the fit.
3. Refit to the cleaned dataset.

The quadratic coefficients encode corrections to the pulsar **period**, **period derivative**, and **reference epoch**.

---

## File Structure

```
.
├── Dedispersion.py
├── Dedispersion and data visualisation.ipynb
├── New Dedispersion and data visualisation.ipynb
├── Searching for Pulsars.ipynb
├── Crab Pulsar Timing.ipynb
├── Make Time of Arrivals.ipynb
├── Obs/
│   ├── B0329+54/
│   ├── B0531+21/              # Crab pulsar observations
│   ├── B1933+16/
│   ├── B2020+28/
│   ├── B2021+51/
│   ├── ASCII/                 # Pulse profiles in ASCII format
│   ├── ASCIItemplates/        # Standard templates
│   └── toas/                  # Computed TOA files
```

---

## Dependencies

- `numpy`
- `scipy`
- `matplotlib`
- `astropy`

Install with:
```bash
pip install numpy scipy matplotlib astropy
```
