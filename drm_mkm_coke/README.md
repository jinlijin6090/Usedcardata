# DRM-MKM with explicit coking steps

Extension of the Ni/CeO2 dry-reforming multisite microkinetic model of
Naveen et al., *ACS Catal.* 2026, 16, 2535 (code: github.com/shenoy1099/DRM-MKM)
with three carbon-carbon coupling / coke steps. The published 17-step model is
reproduced exactly when the coke steps are switched off.

```
mkmc/                 engine (copied from the paper's `mkm` package and extended)
  network.py          20 reactions / 16 species; `use_variant("blocking")` gives 21 / 17
  coke.py             CokeParams + extend_rate_constants + load_model
MKM_parameters.xlsx   the paper's parameter workbook (sheet "Refined params" is used)
coke_tools.py         solve(), drc(), sweep() helpers
s1_barrier_scan.py    sensitivity to Ea_cc, dH_cc, Ea_ox, Ea_dep
s2_condition_sweeps.py pCH4, pCO2, rm, rint, pH2 cofeed, T
s3_drc.py             degree of rate control of DRM rate / coke rate / coke selectivity
s4_sobol.py           Sobol global sensitivity (python s4_sobol.py N PROCS)
s5_blocking.py        scenario B: coke stays on Ni sites, steady state + transient
s6_maps.py            2-D operating maps of coke selectivity
filament.py           scenario C engine: bulk dissolution / diffusion / filament precipitation +
                      surface carbon islands with perimeter-limited growth (site-normalised ODE)
s7_filament.py        scenario C analyses: pure routes, regime map, particle-population model
results/              CSV tables
figures/              PNG figures
```

Run everything (needs numpy scipy pandas matplotlib openpyxl SALib):

```bash
python s1_barrier_scan.py; python s2_condition_sweeps.py; python s3_drc.py
python s4_sobol.py 256 4; python s5_blocking.py; python s6_maps.py; python s7_filament.py
```

Change the coke parameters in `mkmc/coke.py` (`CokeParams`) or pass
`CokeParams(Ea_cc=..., dH_cc=..., Ea_ox=..., Ea_dep=..., dH_g=..., Ea_gas=...)`
to `load_model`.
