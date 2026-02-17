# Modernization Plan Summary

**Project**: Laurant's Thesis → Modern Python Weak-Lensing Pipeline  
**Decision Date**: February 17, 2026  
**Status**: ✓ All recommendations approved, ready to implement  

---

## The Four-Section Review

### **Section 1: Architecture Review** ✓ APPROVED

| Issue | Problem | Decision | Impact |
|-------|---------|----------|--------|
| **1. Monolithic Pipeline** | Rigid linear flow; can't swap components or test subsets | Modular pipeline: DataLoader → JacobianComputer → ObservableComputer → Correlator → Fitter | Can test on sky patches; enables future extensions |
| **2. Scattered Cosmology Params** | Hardcoded in 4+ files; single change requires multi-file edits | Centralized YAML config + CosmologyConfig dataclass | Single source of truth; reproducible |
| **3. Implicit Redshift/ID Coupling** | Parallel arrays (ids, zlist) prone to mismatches | RedshiftSlice dataclass + structured catalog | Type-safe; prevents silent data errors |
| **4. No Test Data Pipeline** | No way to test without full 1.7GB FITS files; slow iteration | Synthetic map generator + small real patches (<100MB) | 10× faster iteration, testable without massive files |

**Cost**: 2-3 weeks | **Benefit**: 10× improvement in testability & iteration speed

---

### **Section 2: Code Quality Review** ✓ APPROVED

| Issue | Problem | Decision | Impact |
|-------|---------|----------|--------|
| **1. Massive DRY Violation** | NFW functions duplicated in 3 files; bug fixes scattered | Unified `lensing_models/nfw_profile.py` with all consolidate functions | Single point of maintenance; easy to add new profiles later |
| **2. Missing Input Validation** | Silent failures on wrong array shapes; bugs hard to trace | Explicit contracts: type checks + shape validation + NaN/Inf detection | Catches errors early; better debugging |
| **3. Mixing Concerns** | Physics + I/O + plotting tangled together; can't test physics alone | Strict separation: Physics → I/O wrapper → Plotting layer (phased refactor) | Each layer independently testable; easy to swap backends |
| **4. Zero Documentation** | Formulas not cited; future students can't verify | LaTeX docstrings + paper citations for every formula | Thesis credibility; science reproducibility |

**Cost**: 4-5 weeks | **Benefit**: Maintainability + scientific rigor + extensibility

---

### **Section 3: Testing Architecture** ✓ APPROVED

| Issue | Solution | Test Count | Benefit |
|-------|----------|-----------|---------|
| **1. No Unit Tests** | pytest suite w/ synthetic fixtures | 60+ (Phase 1) | Regression detection; CI/CD enabled |
| **2. No Reference Values** | Julia→Python HDF5 baseline comparison | Auto-regression | Catch algorithmic changes immediately |
| **3. No Small-Patch Data** | Synthetic generator + real 64×64 corners | 2 fixture types | Millisecond tests; laptop-testable |
| **4. No Perf Baseline** | Timing SLAs (Jacobian <100ms, pipeline <10s) | ~5-10 benchmarks | Auto-detect performance regressions |

**Cost**: 3-4 weeks | **Benefit**: 100+ tests; fast feedback loops; regression safety

**Timeline**: Phase 1 (2 wks), Phase 2 (1 wk), Phase 6 (1 wk)

---

### **Section 4: Modern Libraries & Performance** ✓ APPROVED

| Stage | Current | Python Stack | Rationale |
|-------|---------|--------------|-----------|
| **Jacobian Computation** | Julia @simd loops (hand-optimized) | Numba JIT (speed ≈ Julia) + scipy.ndimage fallback | Preserve performance; clean code |
| **Cosmology Calculations** | Custom code + external C libs | **astropy.cosmology** (industry standard) | Don't reinvent; relies on experts |
| **NFW Lensing Model** | Scattered (3+ files) | Custom NFWProfile class (from Code Quality #1) | Science-specific; single source |
| **2-Point Correlations** | Custom `comp_gs_corr()` (hidden) | **treecorr** (KD-tree, optimized) | 100× faster; handles cosmic variance |
| **I/O Strategy** | Fortran binary + FITS + JLD2 | **HDF5** (single file, compressed, portable) | National standard; compression ~80% |

**Performance Target**: 64×64 Jacobian <100ms, full pipeline <10s

**Cost**: 2-3 weeks | **Benefit**: Industry-standard tools; 10× shorter iteration; eliminates custom C libs

---

## Key Decisions Made

### Architecture
- ✓ Modular pipeline with DataLoader, Physics, I/O, Plotting as separate layers
- ✓ Config-driven (YAML) instead of hardcoded parameters
- ✓ Centralized redshift catalog (dataclass-based)

### Code Quality
- ✓ All NFW/cosmology consolidated to `lensing_models/`
- ✓ Every function has input validation + shaped contracts
- ✓ Full LaTeX docstrings with paper citations
- ✓ Phased separation of I/O/physics/plotting (Phase 1 just consolidation, Phase 2-3 separation)

### Testing
- ✓ 100+ tests across unit/integration/performance tiers
- ✓ Julia baseline comparison for regression detection
- ✓ Synthetic + real small-patch fixtures for fast iteration

### Libraries
- ✓ Numba JIT for Jacobian (matches Julia speed)
- ✓ astropy.cosmology (no reinvention)
- ✓ treecorr for correlations (standard in weak lensing)
- ✓ HDF5 for all I/O (one-time migration from Fortran)

---

## Roadmap at a Glance

```
Phase 0: Setup (Days 1-2)
  └─ Project structure, pytest, CI/CD

Phase 1: Core Physics (Weeks 1-2) ← CRITICAL PATH
  ├─ Cosmology config + astropy wrapper
  ├─ NFW profile consolidation
  ├─ Synthetic + Julia reference data
  ├─ Jacobian (Numba JIT)
  └─ Observable extraction (κ, γ, F, G)
  Status: 60+ unit tests, regression ready

Phase 2: Data & Pipeline (Weeks 3-4)
  ├─ I/O (Fortran binary → HDF5)
  ├─ Data loader
  ├─ Correlator (treecorr)
  └─ Fitter (NFW parameter fitting)
  Status: End-to-end on small data

Phase 3: Orchestration (Week 5)
  ├─ Main entry point
  ├─ Config system
  └─ Integration tests
  Status: Full pipeline runs

Phase 4: Plotting (Week 6)
  ├─ Correlation visualization
  ├─ Map plots
  └─ Publication-quality figures
  Status: Thesis figures ready

Phase 5: Documentation (Week 7)
  ├─ Physics formulas + citations
  ├─ API reference
  ├─ Quickstart guide
  └─ README
  Status: Thesis-ready docs

Phase 6: Validation (Week 8)
  ├─ Julia regression tests
  ├─ Performance profiling
  └─ Final optimization
  Status: Production-ready
```

**Estimated Duration**: 8 weeks (4-5 weeks with 2-3 person team)

---

## Success Criteria

### By End of Phase 1 (Week 2)
- [ ] 60+ unit tests passing
- [ ] Jacobian & observables match Julia reference (<1e-5 relative error)
- [ ] Code coverage >85%

### By End of Phase 3 (Week 5)
- [ ] Full pipeline runs end-to-end on synthetic data
- [ ] Produces: jacobians, correlations, NFW fits
- [ ] Performance: 64×64 map in <10 seconds

### By End of Phase 5 (Week 7)
- [ ] 100+ tests, all passing
- [ ] All formulas documented with citations
- [ ] `pip install -e . && python scripts/run_pipeline.py` works

### By End of Phase 6 (Week 8)
- [ ] Python results match Julia within ~0.1%
- [ ] Full pipeline on real data (10 z-slices)
- [ ] Publication-ready figures + documentation
- [ ] Ready for thesis submission + open-source release

---

## File Structure (Target)

```
lensing-pipeline/
├── lensing/
│   ├── config.py              # YAML config + dataclasses
│   ├── cosmology.py           # astropy wrapper
│   ├── lensing_models/
│   │   ├── nfw_profile.py     # All NFW functions
│   │   ├── jacobian.py        # Numba JIT computation
│   │   └── observables.py     # κ, γ, F, G extraction
│   ├── pipeline/
│   │   ├── data_loader.py
│   │   ├── correlator.py      # treecorr wrapper
│   │   ├── fitter.py
│   │   └── __init__.py        # Orchestrator
│   ├── io/
│   │   ├── fortran_binary.py
│   │   ├── hdf5_backend.py
│   │   └── converters.py
│   └── plotting/
│       ├── correlation_plots.py
│       └── map_plots.py
├── tests/
│   ├── conftest.py
│   ├── fixtures/
│   │   ├── synth_maps.py      # Synthetic generators
│   │   └── reference_outputs.h5  # Julia baseline
│   ├── unit/                  # 60+ unit tests
│   ├── integration/           # 10-15 integration tests
│   └── performance/           # Timing SLAs
├── scripts/
│   ├── run_pipeline.py        # Main entry point
│   ├── convert_fortran_to_hdf5.py
│   ├── validate_against_julia.py
│   └── plot_correlations.py
├── data/
│   ├── config/
│   │   └── thesis.yaml        # Master config
│   └── processed/             # HDF5 outputs
├── docs/
│   ├── physics.md             # Formula compendium
│   ├── api.md                 # Function reference
│   ├── architecture.md
│   └── quickstart.md
├── pyproject.toml
├── requirements.txt
├── MODERNIZATION_ROADMAP.md   # This plan (expanded)
└── README.md
```

---

## Next Steps (Action Items)

### Week 1: Kickoff
- [ ] Create project structure (`pyproject.toml`, `pytest.ini`, directory tree)
- [ ] Set up git, GitHub Actions CI/CD
- [ ] Set up pre-commit hooks (black, mypy, isort)

### Week 1-2: Core Physics (Parallel tasks)
- [ ] **Person A**: Phase 1a-1e (Cosmology → Observables)
  - Config + cosmology
  - NFW consolidation
  - Jacobian (Numba)
  - Observable extraction
  
- [ ] **Person B**: Phase 1c (Synthetic data)
  - Synthetic map generators
  - Julia reference outputs
  - Test fixtures

### Week 2 Checkpoint
```bash
pytest tests/unit/ -v  # Should show ✓✓✓ (60+ tests)
```

### Weeks 3-4: Data Integration
- [ ] **Person B**: Data I/O (Fortran → HDF5) + loader
- [ ] **Person C**: Correlator (treecorr) + Fitter

### Week 5: Orchestration
- [ ] **Person A**: Main script + pipeline config
- [ ] Full end-to-end on synthetic data

---

## Support Resources

**Key Papers** (for docstrings/documentation):
- Bartelmann & Schneider (2001) - Weak lensing theory bible
- Navarro, Frenk & White (1997) - NFW density profile
- Hogg (1999) - Cosmological distance measure conventions
- Schneider et al. (1998) - Lensing observables (κ, γ, flexion)

**Libraries**:
- [astropy.cosmology](https://docs.astropy.org/en/stable/cosmology/)
- [treecorr](https://rmjarvis.github.io/TreeCorr/)
- [Numba JIT](https://numba.readthedocs.io/)
- [h5py](https://docs.h5py.org/)

**Example Real Pipelines**:
- [Describer (LSST)](https://github.com/LSSTDESC/): Weak lensing at scale
- [TreeCorr Examples](https://rmjarvis.github.io/TreeCorr/): 2-point correlation usage

---

## FAQ

**Q: Can we parallelize this across a team?**  
A: Yes! Phase 1a/1c can run in parallel. Phase 2a/2b/2c can overlap. With 3 people, you reduce 8 weeks to ~4 weeks elapsed time.

**Q: What if Numba JIT doesn't match Julia speed?**  
A: Fallback to scipy.ndimage.convolve1d (still >100× faster than naive Python loops). If critical, use Cython (more effort but guaranteed C-level speed).

**Q: Do we need to keep Laurant's Julia code?**  
A: Yes, until Phase 6 validation is complete (need it for regression baseline). After that, archive it.

**Q: Can real data be tested before Phase 3?**  
A: Yes! Use small real patches (extracted in Phase 1c) for integration tests starting Phase 2b.

**Q: What about GPU acceleration?**  
A: Out of scope for now. JAX + GPU could 10× Jacobian computation, but only after Python version is stable (Phase 6+).

**Q: How do we handle thesis deadlines?**  
A: The roadmap is a **minimum viable product**. After Phase 3, you can use the pipeline for thesis work while Phases 4-5 polish documentation.

---

## Sign-Off

**All recommendations approved**: ✓

**Architecture Issues (4)**: All approved  
**Code Quality Issues (4)**: All approved  
**Testing Architecture (4)**: All approved  
**Modern Libraries (4)**: All approved  

**Ready to begin Phase 0 (Setup) on**: Monday, Feb 24, 2026  

---

*For detailed phase breakdown, implementation code outlines, and test examples, see [MODERNIZATION_ROADMAP.md](./MODERNIZATION_ROADMAP.md).*

