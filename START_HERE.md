# 🎯 Modernization Plan: Start Here

**Project**: Laurant's Weak-Lensing Codebase → Modern Python Pipeline  
**Decision Date**: February 17, 2026  
**Status**: ✅ All sections reviewed and approved. Ready to implement.

---

## Which Document Should I Read?

Choose based on your role and time:

### 📋 **If you have 10 minutes** → Read [PLAN_SUMMARY.md](./PLAN_SUMMARY.md)
Quick overview of all decisions, key recommendations, and 8-week timeline.

**What you'll get**: High-level picture, "why" behind each recommendation, success criteria.

---

### 🗺️ **If you have 1-2 hours** → Read [MODERNIZATION_ROADMAP.md](./MODERNIZATION_ROADMAP.md)
Detailed phase-by-phase implementation guide with code outlines, tests, and deliverables.

**What you'll get**: Everything - architecture, dependencies, file structure, critical path, risks.

**Best for**: Project leads, architects, anyone planning resource allocation.

---

### ✅ **If you're implementing** → Use [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md)
Task-by-task checklist with owner assignments, deliverables, and completion criteria.

**What you'll get**: Concrete checkboxes, test requirements, checkpoints, tracking sheet.

**Best for**: Developers executing the work, keeping progress synchronized.

---

## The Plan at a Glance

### What We're Fixing

| Problem | Solution | Benefit |
|---------|----------|---------|
| Scattered Julia code | Modular Python pipeline | 10× easier to test and extend |
| DRY violations | Consolidated modules (NFW, cosmology, etc.) | Single source of truth |
| No tests | 100+ tests (unit + integration) | Regression safety, CI/CD enabled |
| Implicit parameters | Config-driven (YAML) | Reproducible, auditable |
| Custom code | Standard libraries (astropy, treecorr, Numba) | Leverage experts, reduce bugs |

### Timeline

```
Week 1-2:  Core Physics (Cosmology, NFW, Jacobian, Observables) → 60+ tests
Week 3-4:  Data I/O & Pipeline (HDF5, loader, correlator, fitter)
Week 5:    Orchestration (main script, end-to-end on synthetic data)
Week 6:    Plotting (publication-quality figures)
Week 7:    Documentation (formulas, API, quickstart)
Week 8:    Validation (Julia regression, performance tuning)

Total: 8 weeks (or ~4-5 weeks with a 2-3 person team)
```

### Key Decisions ✓ All Approved

**Architecture**:
- Modular pipeline (DataLoader → Physics → I/O → Plotting)
- Config-driven (centralized YAML parameters)
- Structured redshift catalog (dataclass-based)

**Code Quality**:
- Unified NFW module (DRY fix)
- Input validation on all functions
- LaTeX docstrings with paper citations
- Phased separation of concerns

**Testing**:
- 100+ tests (unit + integration)
- Synthetic + real fixtures
- Julia regression baseline
- Performance SLAs (Jacobian <100ms)

**Dependencies**:
- **Numba JIT** for Jacobian (speed ≈ Julia)
- **astropy.cosmology** (industry standard)
- **treecorr** for correlations (weak-lensing specialist library)
- **HDF5** for uniform I/O (compressed, portable)

---

## Getting Started (Next 3 Days)

### Day 1: Understand the Plan
1. Read [PLAN_SUMMARY.md](./PLAN_SUMMARY.md) (15 min)
2. Skim [MODERNIZATION_ROADMAP.md](./MODERNIZATION_ROADMAP.md) phases 1a-1e (30 min)
3. Review this directory structure in your head

### Day 2: Set Up Infrastructure (Phase 0)
- [ ] Create `pyproject.toml`, `pytest.ini`, GitHub Actions workflow
- [ ] Initialize `/lensing`, `/tests`, `/scripts`, `/data` directories
- [ ] First commit to git
- [ ] Verify: `pytest tests/ -v` runs (0 tests, that's fine)

### Day 3: Assign Roles & Schedule
- [ ] Assign Persons A, B, C to phases (see ROADMAP)
- [ ] Schedule kickoff meeting
- [ ] Create GitHub Project board with Phase 1 tasks
- [ ] Kick off Phase 1a (cosmology module)

---

## What Success Looks Like

### After Phase 1 (2 weeks)
```bash
pytest tests/unit/ -v
# → PASSED (60+ tests)
# ✓ All core physics functions implemented
# ✓ Regression tested against Julia baseline
# ✓ Ready for integration testing
```

### After Phase 3 (5 weeks)
```bash
python scripts/run_pipeline.py data/config/thesis.yaml --z-ids 50
# → Complete in <10 seconds
# → Outputs: jacobians/, correlations/, fits/
# ✓ Full pipeline runs end-to-end
# ✓ Ready to start thesis analysis
```

### After Phase 6 (8 weeks)
```bash
pytest tests/ -v --tb=short
# → PASSED (100+ tests)
# Coverage: >85%
# Regression: Python ≈ Julia (<0.1% error)
# Performance: All SLAs met
# ✓ Production-ready, thesis-ready, publication-ready
```

---

## FAQ Quick Answers

**Q: Can we use this for thesis work before Phase 8 is done?**  
A: Absolutely! After Phase 3 (week 5), the pipeline works end-to-end. Phases 4-6 are polish. You can start analyses in week 5 while documentation/validation happens.

**Q: What if we don't have 3 people?**  
A: Timeline stretches but stays manageable:
- 1 person: ~16 weeks (do phases sequentially)
- 2 people: ~8-10 weeks (parallelize prep work)
- 3+ people: ~4-5 weeks (full parallelization)

**Q: Do we keep Julia code?**  
A: Yes, until Phase 6 validation is complete (need it for reference outputs). After that, archive it or keep as documentation.

**Q: Will we use Numba or scipy.ndimage for Jacobian?**  
A: Start with Numba (try to match Julia speed). If Numba disappoints in profiling, fallback to scipy.ndimage.convolve1d (still 100× faster than naive Python loops).

**Q: How do we validate against Julia?**  
A: Pre-generate reference outputs on synthetic maps (Phase 1c), then compare Python results in unit tests (Phase 1d-1e). Later, do full integration tests on real data (Phase 6).

---

## Key Files to Read First

1. **[PLAN_SUMMARY.md](./PLAN_SUMMARY.md)** — 10 minutes, executive summary
2. **[MODERNIZATION_ROADMAP.md](./MODERNIZATION_ROADMAP.md)** — 1 hour, for detailed planning
3. **[IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md)** — Reference while working

Additional reference (already provided):
- `.github/copilot-instructions.md` — Your original design preferences

---

## Contact & Escalation

**Decision already made?** → Proceed with Phase 0 kickoff  
**Questions on strategy?** → Review relevant section in ROADMAP  
**Blocked during implementation?** → Check RISKS section in ROADMAP, consult CHECKLIST  
**Performance issue?** → See Phase 6 (Validation) for profiling guidance  

---

## Success Metrics (Read This Last)

You'll know you've succeeded when:

✅ Code is testable (100+ tests, >85% coverage)  
✅ Science is reproducible (all formulas cited, config-driven)  
✅ Performance is acceptable (Jacobian <100ms, pipeline <10s for small data)  
✅ Pipeline is extensible (easy to add new lensing profiles, new observables)  
✅ Thesis can be written & tested in real-time (no month-long recomputions)  
✅ Code can be published (clean, documented, no technical debt)  

---

**Ready to start? Begin with Phase 0 in [IMPLEMENTATION_CHECKLIST.md](./IMPLEMENTATION_CHECKLIST.md).**

Good luck! 🚀

---

**Plan created**: February 17, 2026  
**Approval status**: ✅ All sections approved  
**Next milestone**: Phase 0 complete by end of week 1

