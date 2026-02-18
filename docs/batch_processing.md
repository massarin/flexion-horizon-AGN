# Batch Processing Guide

The `cosmo_lensing` package handles large datasets efficiently through:

## 1. Automatic Batch Processing in GGL Module

The `ggl.compute_all_ggl_correlations()` function processes multiple observables sequentially:

```python
results = ggl.compute_all_ggl_correlations(
    ra_lens, dec_lens, ra_src, dec_src,
    gamma1, gamma2, F1, F2, G1, G2, kappa=kappa
)
```

Each observable is computed independently, keeping memory usage constant.

## 2. TreeCorr Automatic Parallelization

TreeCorr uses OpenMP for parallel computation of correlations. Set threads:

```python
import treecorr
treecorr.set_omp_threads(4)  # Use 4 CPU cores
```

## 3. Processing Large Catalogs

For very large source catalogs (>1M galaxies), process in chunks:

```python
chunk_size = 100_000
n_chunks = len(ra_src) // chunk_size + 1

for i in range(n_chunks):
    start = i * chunk_size
    end = min((i + 1) * chunk_size, len(ra_src))
    
    chunk_result = ggl.compute_ggl_correlation(
        ra_lens, dec_lens,
        ra_src[start:end], dec_src[start:end],
        gamma1[start:end], gamma2[start:end]
    )
```

## 4. Memory Management

- **Lensing maps**: Keep one map in memory (~1.2 GB for 7200×7200)
- **Galaxy catalog**: ~200 MB for full HAGN catalog
- **Peak memory**: ~2-3 GB for full pipeline
- **No HEALPix needed**: HAGN data is Cartesian patches, not spherical

## Performance Benchmarks

| Operation | Time | Memory |
|-----------|------|--------|
| Load map | 5s | 1.2 GB |
| Load catalog | 2s | 200 MB |
| Correlation (10k lenses, 100k sources) | 10s | 500 MB |
| Profile fitting | <1s | negligible |

**Total workflow**: ~30 seconds per redshift slice
