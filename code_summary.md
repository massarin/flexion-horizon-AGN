# Code Summary: Weak Gravitational Flexion in Horizon-AGN Simulation

## Overview

This codebase implements the analysis pipeline for measuring weak gravitational lensing signals (convergence κ, shear γ, and flexion F/G) from the Horizon-AGN hydrodynamical simulation. The thesis measures galaxy-galaxy lensing correlations at various redshifts and compares signals between central and satellite galaxies, fitting results to SIS and NFW mass profile models.

**Primary Language**: Julia  
**Data Source**: Horizon-AGN simulation lightcone with raytracing  
**Main Output**: Correlation functions for 7 lensing fields (κ, γ₁, γ₂, F₁, F₂, G₁, G₂)

---

## Core Data Flow

```
Input Data
    ├── Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits  → Galaxy catalog (positions, masses, redshifts)
    └── Binary deflection maps (α₁, α₂)             → Ray-traced deflection fields
         ↓
Lensing Map Generation (rebuild_*.jl)
    ├── Read deflection maps → Compute Jacobian (alpha2jac)
    └── Generate lensing maps: κ, γ₁, γ₂, F₁, F₂, G₁, G₂
         ↓
Correlation Computation (correlation_function.jl)
    ├── Select galaxies by: redshift, mass, central/satellite
    ├── Compute radial averages around each galaxy
    └── Apply Jackknife error estimation
         ↓
Model Fitting &amp; Visualization (correlation_plots.jl)
    ├── Fit NFW and SIS profiles
    ├── Extract halo mass and concentration
    └── Generate plots for thesis
```

---

## File-by-File Breakdown

### Core Infrastructure

#### `core.jl` (6 lines)
- **Purpose**: Interface to external softlens C library
- **Key Function**: Defines `sllib` path to libsoftlens shared library
- **Status**: Minimal wrapper; actual cosmology functions in cosmo.jl

#### `cosmo.jl` (131 lines)
- **Purpose**: Cosmological calculations using softlens library + Julia Cosmology.jl
- **Key Functions**:
  - `set_cosmology()`: Configure ΛCDM parameters (H₀=70.4, Ωₘ=0.3, ΩΛ=0.7)
  - `dda()`, `dlsdos()`: Angular diameter distances and ratios
  - `arcsec2kpc()`: Physical scale conversion
  - `Δ_vir()`, `ρ_crit()`: Virial overdensity and critical density
- **Used By**: All analysis scripts
- **Note**: Defines flat ΛCDM matching Horizon-AGN parameters

#### `nfw_tk.jl` (104 lines)
- **Purpose**: NFW profile utilities and mass-concentration relations
- **Key Functions**:
  - `mc2ksrs()`: Convert (mass, concentration) → (κₛ, rₛ) NFW parameters
  - `ksrs2mc()`: Inverse conversion (κₛ, rₛ) → (mass, concentration)
  - `c_m_K16()`: Concentration-mass relation from Klypin et al. 2016
- **Used By**: correlation_function.jl for converting fit parameters to physical quantities
- **Note**: Handles both scalar and vectorized conversions

### Data Processing Core

#### `raytrace_tk.jl` (1734 lines) ⭐ **CENTRAL FILE**
- **Purpose**: Complete toolkit for raytracing, Jacobian computation, and correlation analysis
- **Key Components**:

**1. I/O Functions** (lines 67-160)
  - `read_bin()`: Read Fortran binary deflection maps (chunked for memory efficiency)
  - `write_bin()`: Write deflection fields
  - `average_law()`: Binned averaging utility

**2. Jacobian Computation** (lines 218-530)
  - `alpha2jac()`: **Core function** - computes all derivatives of deflection field
    - Input: α(x,y) deflection map [nx × ny × 2]
    - Output: Jacobian matrix [nx × ny × 12] containing:
      - ∂₁α₁, ∂₁α₂, ∂₂α₁, ∂₂α₂ (first derivatives)
      - ∂₁₁α₁, ∂₁₁α₂, ∂₂₂α₁, ∂₂₂α₂ (second derivatives)
      - ∂₁₂α₁, ∂₁₂α₂, ∂₂₁α₁, ∂₂₁α₂ (cross derivatives for flexion)
    - Uses finite differences with periodic boundaries
  - Alternative versions: `alpha2jac_old()`, `alpha2jac_norot_commder()`, `alpha2jac_flexion_w_cross_terms()`
  
**3. Lensing Map Extraction** (lines 648-758)
  - `jac2kappa()`: κ = ½(∂₁α₁ + ∂₂α₂)
  - `jac2gamma1()`: γ₁ = ½(∂₁α₁ - ∂₂α₂)
  - `jac2gamma2()`: γ₂ = ½(∂₁α₂ + ∂₂α₁)
  - `jac2F1()`: F₁ = ½(∂₁₁α₁ + ∂₁₂α₂)
  - `jac2F2()`: F₂ = ½(∂₁₁α₂ - ∂₁₂α₁)
  - `jac2G1()`: G₁ = ½(∂₁₁α₁ - ∂₁₂α₂)
  - `jac2G2()`: G₂ = ½(∂₁₁α₂ + ∂₁₂α₁)
  - `jac2rot()`: Rotational component check (should be zero)
  - `jac2muinv()`: Inverse magnification

**4. Galaxy-Galaxy Lensing** (lines 1657-1734)
  - `comp_gs_corr()`: **Main correlation function**
    - Reads deflection map for given lens plane
    - Computes all 7 lensing fields
    - For each galaxy in catalog:
      - Draws concentric radial bins
      - Averages tangential components (κ, γ₊, F₊, G₊) in each bin
    - Returns: Sres[5, nr, ng] array (5 quantities × radial bins × galaxies)
  - `winc()`: Window function for radial binning
  - `get_bbox2()`: Compute bounding box around galaxy

**5. Utility Functions**
  - `rebin()`: Downsample maps by averaging
  - `deriv!()`: Finite difference derivatives
  - `proj_rt()`: Tangential/cross projection of lensing quantities

**Status**: Production code, heavily optimized with @threads, @simd

---

### Analysis Scripts

#### `correlation_function.jl` (419 lines) ⭐ **MAIN ANALYSIS**
- **Purpose**: End-to-end correlation measurement with error estimation
- **Workflow**:

**1. Configuration** (lines 1-142)
  - Redshift list: z = [0.21, 0.35, ..., 3.91] (10 lens planes)
  - Plane IDs: [50, 100, 150, ..., 450]
  - Mass selection: 0=low, 1=intermediate, 2=high (top 500 galaxies)
  - Galaxy type: 0=central, 1=satellite
  - Error method: 1=jackknife, 0=bootstrap

**2. Galaxy Selection** (lines 148-240)
  - Load FITS catalog
  - Redshift cut: z_lens/2 ± 0.1 (optimal lensing efficiency)
  - Central/satellite split: IDH_SUB == IDH_MAIN → central
  - Mass ranking: selects top/middle/bottom 500 galaxies by MTOTH_MAIN

**3. Correlation Computation** (line 243)
  - Calls `comp_gs_corr()` with selected galaxy positions
  - Radial bins: rmax=50 arcmin, nr=50 bins
  - Saves results: Data/corr_data/Sres_rbin_{spec}_id_{id}_arcm50.jld2

**4. Error Estimation** (lines 295-346)
  - **Jackknife**: 5000 iterations, exclude 1% of galaxies each time
  - Computes variance: σ² = (n-1) × mean((x_i - mean)²)
  - Generates error bars for all 4 lensing components

**5. Model Fitting** (lines 392-418)
  - Fits NFW convergence profile: κ(r) = 2κₛrₛ²(1-F(r/rₛ))/(r²-rₛ²)
  - Uses LsqFit.jl curve_fit
  - Converts (κₛ, rₛ) → (mass, concentration) via ksrs2mc()
  - Initial guess: p0 = [1.0, 0.01]

**6. Visualization** (lines 365-390)
  - Generates scatter plots with error bars
  - Separate plots per component: κ, γ, F, G
  - Log-log and semi-log variants
  - Saves to corr_plots/{component}/

**Reported in Thesis**: 
  - Figures 4.4, 4.5 (correlation plots for all components)
  - Figures 4.6, 4.7 (NFW fits for central/satellite galaxies)

#### `correlation_plots.jl` (467 lines)
- **Purpose**: Plotting-only version (assumes correlations already computed)
- **Difference from correlation_function.jl**: Skips comp_gs_corr(), loads JLD2 files
- **Use Case**: Regenerate plots after tweaking visualization parameters
- **Status**: Actively used for thesis figures

#### `fitting_correlation.jl` (115 lines)
- **Purpose**: Standalone fitting experiments
- **Content**: Tests various NFW/SIS fitting approaches
- **Status**: Experimental/exploratory, not used in final thesis

---

### Map Generation Scripts

These scripts generate FITS files containing individual lensing maps for visualization. They read deflection maps, compute derivatives, and save rebinned results.

#### `rebuild_maps.jl` (105 lines)
- **Purpose**: Generate full lensing derivative maps for all planes
- **Output**: HAGN-lensingderiv2ord_{id:04d}.fits containing:
  - kappa, gamma1, gamma2, F1, F2, G1, G2
  - rot (rotation check), crossderivatives (a1, a2)
  - imu (inverse magnification)
- **Rebinning**: rf=5 (reduces 36000² → 7200² for storage)
- **Used For**: Figure 4.1, 4.2, 4.3 (map visualizations in thesis)

#### `rebuild_{component}_maps.jl` (61-95 lines each)
- **Files**: rebuild_convergence_maps.jl, rebuild_gamma1_maps.jl, rebuild_gamma2_maps.jl,
  rebuild_F1_maps.jl, rebuild_F2_maps.jl, rebuild_G1_maps.jl, rebuild_G2_maps.jl,
  rebuild_shear_maps.jl, rebuild_flexion_maps.jl, rebuild_rot_maps.jl
- **Purpose**: Generate individual component FITS files
- **Pattern**: 
  1. Read deflection map
  2. Compute Jacobian with alpha2jac()
  3. Extract specific component (e.g., jac2kappa())
  4. Rebin by factor rf=5
  5. Save to lensing_maps/{component}/
- **Status**: Utility scripts for creating publication-quality maps
- **Not in Thesis**: Individual map files not directly used, only full maps from rebuild_maps.jl

#### `rebuild_map_differences.jl` (80 lines)
- **Purpose**: Compute differences between different Jacobian implementations
- **Use Case**: Validation - checking consistency between alpha2jac variants
- **Status**: Testing/verification only

---

### Model Comparison Scripts

#### `correlation_plots_sis_model.jl` (690 lines)
- **Purpose**: Fit Singular Isothermal Sphere (SIS) profiles
- **SIS Model**: κ(r) ∝ 1/r (simpler than NFW)
- **Content**: Similar structure to correlation_function.jl but with SIS fitting
- **Status**: Exploratory - NFW fits preferred for thesis (better match to simulation)

#### `correlation_plots_sis_model_2.jl` (682 lines)
- **Purpose**: Variant of SIS fitting with different parameters
- **Status**: Iteration/experimentation

---

### Testing &amp; Utilities

#### `test_some_maps.jl` (124 lines)
- **Purpose**: Quick tests for map generation pipeline
- **Content**: Loads subset of planes, generates sample maps
- **Status**: Development utility

#### `juliatest.jl` (296 lines)
- **Purpose**: General Julia syntax and function testing
- **Status**: Scratch file

#### `image_distort.jl` (201 lines)
- **Purpose**: Apply lensing distortions to synthetic galaxy images
- **Content**: Forward lensing (image → source plane transformation)
- **Status**: Exploration - not used in thesis (thesis uses weak lensing statistics, not image reconstruction)

---

## Data Files

### Input Data

#### `Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits`
- **Format**: FITS binary table
- **Size**: Large catalog (millions of galaxies)
- **Key Columns** (from Table 3.2 in thesis):
  - `RA_IMG`, `DEC_IMG`: Lensed sky positions (degrees)
  - `z_true`: True redshift
  - `MTOTH_MAIN`: Total mass of main halo (solar masses)
  - `MTOTH_SUB`: Total mass of subhalo
  - `IDH_MAIN`, `IDH_SUB`: Halo identifiers (used for central/satellite split)
  - `KAPPA`, `GAMMA1`, `GAMMA2`: Pre-computed lensing (not used - recomputed from deflection)
  - `MAGNIF`: Magnification factor

#### Deflection Maps (binary format)
- **Location**: Referenced by `full_bin_name()` in raytrace_tk.jl
- **Format**: Fortran unformatted binary (read via FortranFiles.jl)
- **Structure**:
  - Header: size [nx, ny], dummy array, source redshift
  - Data: α₁(x,y) and α₂(x,y) deflection components [Float32]
  - Size: 36000 × 36000 pixels per plane
- **Planes**: 10 planes at IDs [50, 100, 150, 200, 250, 300, 350, 400, 417, 450]
  - Redshifts: z = [0.21, 0.35, 0.51, 0.73, 1.02, 1.42, 1.82, 2.46, 2.79, 3.91]

### Output Data

#### `Data/corr_data/Sres_rbin_*.jld2`
- **Format**: JLD2 (Julia serialized arrays)
- **Files**: One per (plane, mass_bin, galaxy_type) combination
- **Naming**: `Sres_rbin_{highmass|intermediatemass|lowmass}_{central|satellite}_alpha2jac_final_id_{id:04d}_arcm50.jld2`
- **Content**:
  - `Sres`: [5 × nr × ng] array
    - Dimension 1: [G₊, F₊, κ, γ₊, pixel_count]
    - Dimension 2: Radial bins (50 logarithmically spaced)
    - Dimension 3: Individual galaxies
  - `rbin`: Bin edges in arcminutes

#### `lensing_maps/` (FITS files)
- **Subdirectories**: all/, convergence_maps/, gamma1_maps/, etc.
- **Format**: FITS image extensions
- **File**: HAGN-lensingderiv2ord_{id:04d}.fits or component-specific
- **Size**: ~1.2 GB per full map file (7200² × 10 extensions)

#### `corr_plots/` (PNG plots)
- **Structure**:
  - all/: Combined plots (κ, γ, F, G on same axes)
  - kappa/, gamma/, F/, G/: Individual component plots
  - fits/: Model fitting results
- **Naming**: `{component}_{log|nolog}_{satellite|central}_{highmass|intermediatemass|lowmass}_alpha2jac_final_{arcmin}_arcm_id_{id:04d}.png`

---

## What's Actually in the Thesis

### Reported Results

✅ **Used and Reported**:
1. **Lensing Maps** (Chapter 4.2):
   - Figure 4.1: Convergence map (z≈1.016) from rebuild_maps.jl
   - Figure 4.2: Shear fields γ₁, γ₂ showing spin-2 symmetry
   - Figure 4.3: Flexion fields F₁, F₂, G₁, G₂ showing spin-1/spin-3 symmetry

2. **Correlation Functions** (Chapter 4.3):
   - Figure 4.4: κ, γ, F, G correlations vs. radius (arbitrary planes)
   - Figure 4.5: Same correlations with semi-log plotting
   - Measured for all 10 lens planes
   - Separated by mass bins (high/intermediate/low)
   - Computed via correlation_function.jl → comp_gs_corr()

3. **Central vs. Satellite** (Chapter 4.4):
   - Comparison of 1-halo and 2-halo term behavior
   - Satellite galaxies show softer slope, "bump" at intermediate scales
   - Computed via correlation_function.jl with satellite=0 vs satellite=1

4. **Model Fitting** (Chapter 4.5):
   - Figure 4.6: NFW fits for central galaxies (κ, γ, F, G)
   - Figure 4.7: NFW fits for satellite galaxies
   - Fitted with LsqFit.jl in correlation_function.jl
   - Cutoff radius: ~800 kpc (1-halo term only)
   - Extracted (mass, concentration) parameters via ksrs2mc()

5. **Theoretical Framework** (Chapter 2):
   - Equations for κ, γ, F, G as derivatives of deflection field
   - Tangential projections (Eq. 3.8)
   - NFW profile equations (Eq. 3.20-3.25)

6. **Code Listings** (Appendix A, B):
   - Appendix A: alpha2jac() function (lines 2672-2800+ in thesis)
   - Appendix B: Correlation computation workflow (referenced but not fully shown)

### Partially Discussed

⚠️ **Mentioned but Limited Analysis**:
- Jackknife error estimation (Chapter 4.1.2): Methodology explained, implemented in correlation_function.jl
- Julia performance rationale (Chapter 4.1.1): Why Julia was chosen
- SIS model: Mentioned but not shown in figures (correlation_plots_sis_model.jl exists but not used)

### Not Reported

❌ **Computed but Not in Thesis**:
1. **Cross-components** (κ×, γ×, F×, G×):
   - Computed in comp_gs_corr() but not plotted
   - Should be zero (test for systematics) - computed in test_maps/rot/ etc.

2. **Individual component maps**:
   - rebuild_{gamma1,gamma2,F1,F2,G1,G2}_maps.jl generate separate FITS
   - Only composite map shown (Figure 4.1-4.3)

3. **Alternative Jacobian methods**:
   - alpha2jac_norot_commder(), alpha2jac_flexion_w_cross_terms() exist
   - rebuild_map_differences.jl compares them
   - No discussion of differences in thesis

4. **Bootstrap error estimates**:
   - Implemented (error_estimator=0) but Jackknife preferred

5. **Low/intermediate mass bins**:
   - Computed but only high-mass results shown in figures

6. **Chapter 5 (Euclid/SourceExtractor++/ONNX)**:
   - Discussed as "future work" but no code in this repository
   - Intended to apply flexion measurements to real Euclid data

7. **Image distortion** (image_distort.jl):
   - Not used in thesis at all

---

## Technical Details

### Coordinate Systems
- **Sky coordinates**: RA/DEC in degrees
- **Map coordinates**: Cartesian pixels with flat-sky approximation
- **Pixel scale**: dx = 1.0 arcmin (configurable)
- **Conversion**: 
  - Angular → comoving distance via Cosmology.jl
  - Comoving → physical via scale factor

### Numerical Methods
- **Derivatives**: 5-point finite differences (deriv! function)
- **Rebinning**: Simple averaging over rf×rf blocks
- **Radial bins**: Logarithmically spaced from 0.05 to 50 arcmin
- **Tangential projection**: 
  ```julia
  κ₊ = κ
  γ₊ = -Re(γ × e^(-2iφ))
  F₊ = -Re(F × e^(-iφ))
  G₊ = -Re(G × e^(-3iφ))
  ```

### Performance Optimizations
- **Multithreading**: `Threads.@threads` for per-galaxy loops
- **SIMD**: `@simd` for inner derivative loops
- **Memory management**: Explicit GC.gc(true) after large arrays
- **Chunked I/O**: read_bin() processes large files in 1M element chunks
- **Shared arrays**: Considered but not used (comment in raytrace_tk.jl)

### Dependencies
- **Core**: Julia 1.x (likely 1.6+)
- **Scientific**: Cosmology, Interpolations, Roots, Polynomials, Unitful, UnitfulAstro
- **Data I/O**: FITSIO, FortranFiles, JLD2
- **Plotting**: Plots (likely GR backend)
- **Fitting**: LsqFit
- **Utilities**: Printf, Statistics, DataStructures, ProgressMeter, Random, ForwardDiff
- **External**: libsoftlens (C library, path in ENV["SOFTLENS_DIR"])

---

## Execution Order

### To Reproduce Thesis Results:

1. **Generate lensing maps** (optional, for visualization):
   ```julia
   include("rebuild_maps.jl")
   ```
   - Outputs: lensing_maps/all/HAGN-lensingderiv2ord_*.fits
   - Used for: Figures 4.1, 4.2, 4.3

2. **Compute correlations**:
   ```julia
   # Edit configuration in correlation_function.jl:
   massindic = 2      # 0=low, 1=intermediate, 2=high mass
   satellite = 0      # 0=central, 1=satellite
   error_estimator = 1 # 1=jackknife
   
   include("correlation_function.jl")
   ```
   - Outputs: Data/corr_data/Sres_rbin_*.jld2
   - Used for: All correlation plots and fits

3. **Generate plots** (or use plots from step 2):
   ```julia
   # Same configuration as step 2
   include("correlation_plots.jl")
   ```
   - Outputs: corr_plots/{component}/*.png
   - Used for: Figures 4.4, 4.5, 4.6, 4.7

### Iteration Notes:
- Steps 2-3 need to run for each combination:
  - 2 galaxy types (central, satellite)
  - 3 mass bins (though only high-mass in thesis)
  - = 6 configurations minimum
- Each run takes significant time (hours) due to:
  - Reading 36000² deflection maps
  - Computing 12-component Jacobian
  - Per-galaxy radial averaging (500 galaxies × 50 bins)
  - Jackknife resampling (5000 iterations)

---

## Issues &amp; Observations

### Code Quality

**Strengths**:
- Well-optimized core routines (alpha2jac, comp_gs_corr)
- Proper memory management for large arrays
- Extensive use of Julia performance patterns (@inbounds, @simd)

**Weaknesses**:
- Heavy code duplication across correlation_plots*.jl files
- Hardcoded paths (amalgam vs local machine conditionals)
- Configuration via global variables rather than arguments
- Minimal documentation/comments
- Some French comments mixed with English
- "Crash zone" comment in raytrace_tk.jl (line 215)

### Unclear Elements

1. **Binary deflection map source**: 
   - Where do the original deflection maps come from?
   - Mentioned they're from C. Gouin's raytracing work (reference [6])
   - No generation code in this repository

2. **softlens library**:
   - External C library, minimal documentation
   - Used for cosmology but most functions reimplemented in Julia

3. **Alternative Jacobian functions**:
   - Three variants exist (alpha2jac, alpha2jac_norot_commder, alpha2jac_flexion_w_cross_terms)
   - No clear explanation of which is "correct" or why alternatives exist
   - Production code uses standard alpha2jac()

4. **Fitting cutoff radius**:
   - 800 kpc cutoff mentioned but not clearly justified
   - Separating 1-halo vs 2-halo term is more art than science

5. **Mass values**:
   - Some hardcoded mean masses in correlation_plots.jl (lines 49-52)
   - Not clear if these are used or just reference values

### Missing Validation

- No unit tests
- No comparison with analytical models (except fitting)
- Cross-component (systematic) checks computed but not shown
- No discussion of convergence with spatial resolution

### File Organization

- No clear separation of library vs. scripts
- Backup files in main directory (correlation_plots_sis_model_2.jl)
- Results mixed with code (corr_plots/, lensing_maps/ could be in separate data directory)

---

## Recommendations for Future Work

1. **Refactor for reusability**:
   - Extract configuration to TOML/YAML files
   - Create proper module structure (module LensingAnalysis)
   - Consolidate duplicate code

2. **Add documentation**:
   - Docstrings for all public functions
   - README with setup instructions
   - Example notebooks

3. **Improve reproducibility**:
   - Package environment (Project.toml with explicit versions)
   - Scripts to download/prepare data
   - Makefile or pipeline script

4. **Validation suite**:
   - Unit tests for derivatives
   - Compare with analytical NFW predictions
   - Verify cross-components are consistent with zero

5. **Clean up**:
   - Remove commented code
   - Delete obsolete files (juliatest.jl, backup/)
   - Standardize naming conventions

6. **Performance**:
   - Profile critical sections
   - Consider GPU acceleration for derivative computation
   - Parallelize across multiple lens planes

---

## Summary: What Actually Matters

If you're trying to understand this thesis:

**Focus on**:
1. `raytrace_tk.jl` - lines 218-758 (Jacobian computation, lensing extraction)
2. `correlation_function.jl` - lines 148-243 (galaxy selection, correlation measurement)
3. `correlation_function.jl` - lines 295-418 (error estimation, NFW fitting)

**Ignore**:
- SIS model files (not in thesis)
- Alternative Jacobian functions (not used)
- Image distortion code (different analysis)
- Test files and backup directories

**Key Insight**:
The entire analysis boils down to:
1. Take derivatives of deflection field → get 7 lensing maps
2. Average lensing maps in radial bins around selected galaxies
3. Fit averaged signal to NFW profile → extract mass
4. Compare central vs satellite, different redshifts

The complexity is in the numerical implementation (36000² maps, 500 galaxies, 50 bins, 10 planes) and error estimation (5000 jackknife samples), not in the conceptual framework.
