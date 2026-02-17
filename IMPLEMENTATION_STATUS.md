# Implementation Status Tracker

**Plan Created**: 2026-02-17  
**Status**: Ready for Implementation  
**Full Plan**: See `MODERNIZATION_PLAN.md`

## Quick Summary

Converting scattered Julia thesis code (22 files, ~2500 duplicate lines, no tests) to production-ready Python pipeline with:
- 6 focused modules (io, derivatives, observables, correlations, nfw, cosmology)
- Industry-standard libraries (NumPy, SciPy, Astropy, TreeCorr, HEALPix)
- Comprehensive test suite (unit/integration/validation/production)
- Snakemake workflow for reproducibility
- HPC-ready with full provenance tracking

## Phase Checklist

### Phase 1: Foundation (Week 1-2) - NOT STARTED
- [ ] Create package structure (`cosmo_lensing/`)
- [ ] Set up environment (conda/pip, requirements.txt)
- [ ] Port I/O functions (io.py)
- [ ] Port cosmology (cosmology.py, validate vs libsoftlens)
- [ ] Port derivatives (derivatives.py)
- [ ] Port observables (observables.py with citations)
- [ ] Set up pytest infrastructure

**Deliverable**: Core library can read deflection, compute observables  
**Test**: Process one deflection field, save FITS map

### Phase 2: Synthetic Data & Validation (Week 3) - NOT STARTED
- [ ] Implement synthetic data generators
- [ ] Write 10+ analytic validation tests
- [ ] Implement hierarchical test configs (unit/debug/validation/production)
- [ ] Test on synthetic data at all scales

**Deliverable**: Test suite with all tests passing  
**Test**: `pytest tests/ --config unit` completes in < 10 seconds

### Phase 3: Correlation Functions & TreeCorr (Week 4) - NOT STARTED
- [ ] Implement TangentialCorrelation class
- [ ] Implement NFWProfile class
- [ ] Write TreeCorr validation tests (< 1% agreement required)
- [ ] Port fitting code

**Deliverable**: Correlations match TreeCorr within 1%  
**Test**: Run on synthetic NFW, compare to analytic γ_t(r)

### Phase 4: Production Pipeline (Week 5-6) - NOT STARTED
- [ ] Create Snakefile workflow
- [ ] Create config.yaml template
- [ ] Implement production scripts (CLI interfaces)
- [ ] Add provenance tracking
- [ ] Test workflow (dry run → local → cluster)

**Deliverable**: Fully automated workflow  
**Test**: `snakemake --config scale=validation` completes successfully

### Phase 5: Scaling & Optimization (Week 7) - NOT STARTED
- [ ] Implement HEALPix support
- [ ] Add parallel processing
- [ ] Optimize performance (profile + optimize hotspots)
- [ ] Memory optimization

**Deliverable**: Can process NSIDE=2048 HEALPix (50M pixels)  
**Test**: Production scale run completes in < 1 hour

### Phase 6: Final Validation & Documentation (Week 8) - NOT STARTED
- [ ] Run full production pipeline (10 redshift slices)
- [ ] Compare results to Julia code (< 5% difference)
- [ ] Write comprehensive documentation
- [ ] Code cleanup (delete Julia, linting)
- [ ] Create release (v1.0, Zenodo DOI)

**Deliverable**: Publication-ready code + docs  
**Test**: External user can reproduce thesis results from scratch

## Current Issues (All Approved for Fix)

### Architecture (4 issues)
1. ✅ Monolithic toolkit → 6 focused modules
2. ✅ libsoftlens.so dependency → astropy.cosmology
3. ✅ No spatial decomposition → HEALPix hierarchy
4. ✅ Coupled data/viz → separate layers

### Code Quality (4 issues)
5. ✅ DRY violation (2500 duplicate lines) → single parameterized script
6. ✅ Zero error handling → explicit validation
7. ✅ Undocumented formulas → citations + unit tests
8. ✅ Dead SIS code (50KB) → delete + refactor correlations

### Performance (4 issues)
9. ✅ No test mode → hierarchical test pyramid
10. ✅ No validation vs TreeCorr → cross-validation tests
11. ✅ Research vs production → Snakemake pipeline (approved: production)
12. ✅ Inefficient jackknife → spatial jackknife with TreeCorr

## Key Libraries

| Purpose | Library | Critical? |
|---------|---------|-----------|
| Arrays | NumPy ≥1.24 | Yes |
| Scientific | SciPy ≥1.10 | Yes |
| Astronomy | Astropy ≥5.0 | Yes |
| **Correlations** | **TreeCorr ≥4.3** | **CRITICAL** |
| Spatial | HEALPix ≥1.16 | Yes |
| Workflow | Snakemake ≥7.0 | Yes |
| Testing | pytest ≥7.0 | Yes |

## Next Steps for Dev Agent

1. **Start with Phase 1, Task 1**: Create package structure
   ```bash
   mkdir -p cosmo_lensing tests scripts workflow
   touch cosmo_lensing/__init__.py
   ```

2. **Set up Python environment**
   ```bash
   python -m venv venv
   source venv/bin/activate
   pip install numpy scipy astropy treecorr healpy pytest
   ```

3. **Begin porting io.py**: Start with `read_deflection_field()` function

4. **Reference files**:
   - Original: `raytrace_tk.jl` lines 77-110 (read_bin function)
   - Target: `cosmo_lensing/io.py`
   - Add error handling as specified in Issue #6

## Files to Delete (After Backup)

**Before deleting, create archive:**
```bash
mkdir -p archive/julia_original
cp *.jl archive/julia_original/
git add archive/
git commit -m "Archive original Julia code before Python conversion"
```

**Then delete:**
- All 22 `.jl` files (move to archive)
- `correlation_plots_sis_model.jl` (dead code)
- `correlation_plots_sis_model_2.jl` (dead code)
- 9 duplicate `rebuild_*_maps.jl` scripts

**Keep:**
- `papers/` directory (references)
- `Data/` directory (input data)
- `MODERNIZATION_PLAN.md` (this plan)

## Success Criteria

**Must achieve ALL of these:**
- [ ] All 12 issues resolved (no exceptions)
- [ ] 100% Julia functionality ported
- [ ] Test suite > 80% coverage, all passing
- [ ] TreeCorr validation < 1% difference
- [ ] Fitted masses match Julia < 10% (systematic OK if documented)
- [ ] Pipeline 10× faster (TreeCorr optimization)
- [ ] Memory efficient (32GB laptop for debug scale)
- [ ] Fully reproducible (same input → same output within 0.1%)
- [ ] External user can reproduce thesis from scratch
- [ ] Publication-ready (clean, documented, tested)

## Estimated Timeline

- **Minimum (expert)**: 6 weeks
- **Realistic (part-time)**: 3 months  
- **Conservative (first-time)**: 4 months

**Critical path**: TreeCorr validation (Issue #10) must match within 1%

---

**Status**: ✅ Plan approved, ready for implementation  
**Next**: Dev agent should start Phase 1, Task 1
