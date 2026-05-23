# Literature Benchmark Notes

This folder collects public-reference benchmark data that can be used to check the OPB fatigue workflow at a screening level.

The current public references do not provide full raw time histories of mooring tension, OPB angle, and IPB angle. Because of that, the repository cannot yet reproduce a published damage value from exactly the same time series. The benchmark checks therefore focus on:

- confirming that code constants match the published NI604-style SCF and S-N values;
- confirming that the published fatigue-damage ordering is reproduced as a trend target;
- documenting which extra data would be needed for a strict one-to-one validation.

## Primary Public Reference

Wang et al. (2024), "An Out-of-Plane Bending Fatigue Assessment Approach for Offshore Mooring Chains Considering the Real-Time Updating of Interlink Bending Stiffness", Journal of Marine Science and Engineering.

Link: https://www.mdpi.com/2077-1312/12/1/131

Useful public data:

- Figure 14 describes the same high-level workflow used here: global motion/tension, top-chain interlink angles and moments, hotspot stresses, rainflow counting, S-N damage, and Miner summation.
- Table 6 gives studless-chain hotspot SCFs for TT/IPB/OPB combined stress calculation.
- Table 7 gives S-N parameters, including free corrosion in seawater with `log10(K) = 12.436` and `m = 3`.
- Table 8 gives total fatigue damage and fatigue life for mooring line #2.
- Appendix Table A4 gives case-by-case fatigue damage for mooring line #2.

## Secondary Public References

Bureau Veritas NI604 official page:

https://marine-offshore.bureauveritas.com/ni604-fatigue-top-chain-mooring-lines-due-plane-and-out-plane-bendings

This page identifies NI604 as guidance for combined top-chain fatigue under tension, in-plane bending, and out-of-plane bending.

Engineering Failure Analysis 2019 OPB experiment paper:

https://www.sciencedirect.com/science/article/pii/S1350630718316443

The abstract reports OPB fatigue tests and stress-based fatigue estimates within factor-of-three boundaries. Full test tables may require journal access, so this is currently a qualitative reference only.

## Current Automated Checks

Run:

```powershell
python validation/check_mdpi_trends.py
```

The script checks:

- `STUDLESS_SCF` values used by `ni604_opb_fatigue.py`;
- default S-N values in `SNCurve`;
- Table 8 damage ordering: `Hotspot C > Hotspot B > Hotspot A > TT`;
- Table 8 life ordering and its relation to damage using the published safety factor of 10;
- Appendix Table A4 subset governing-hotspot patterns.

## What Would Be Needed for Strict Reproduction

For strict numerical validation, a public dataset would need:

- tension time series;
- OPB and IPB interlink angle or moment time series;
- exact chain diameter, MBL, pretension, corrosion environment, and safety factors;
- sea-state duration and occurrence probability;
- exact rainflow convention and stress-range/amplitude convention.

With those inputs, `ni604_opb_fatigue.py` can be run directly and compared against the published short-term or long-term damage.
