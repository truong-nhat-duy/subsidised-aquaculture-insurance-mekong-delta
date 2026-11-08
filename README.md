# Why Subsidised Aquaculture Insurance Does Not Reach Farmers

**Climate hazards, trigger design and risk-based pricing in Vietnam's Mekong Delta**

Research code, LaTeX sources and compiled papers for a study of Vietnam's premium-subsidised agricultural insurance for shrimp and pangasius (Decree 58/2018/ND-CP; Decisions 22/2019/QD-TTg, 03/2021/QD-TTg and 13/2022/QD-TTg).

**Website:** https://truong-nhat-duy.github.io/subsidised-aquaculture-insurance-mekong-delta/

## Research question

Premiums are subsidised by group (up to 90% for poor and near-poor households, 20% for others) but are almost flat across hazard zones. Why does cover still barely reach farmers, and would hazard-based trigger design and risk-based pricing do better?

## Evidence summary

- Data: VHLSS household surveys 2002-2024 (restricted; not included), plus climate and hazard layers (ERA5, IBTrACS, CHIRPS/SPEI).
- First stage is close to zero: 0-1% of fish-farming households in pilot provinces report an insurance payout.
- Synthetic difference-in-differences: programme effect 0.016 (SE 0.066), statistically indistinguishable from zero.
- Basis risk of candidate triggers (ROC/AUC): 0.47-0.63, i.e. weak to moderate discrimination.
- Marine heatwave days rise from roughly 7 to roughly 53 per year in the series studied.
- Scope: five provinces (Ben Tre, Tra Vinh, Soc Trang, Bac Lieu, Ca Mau); natural-disaster perils only, no disease cover.

All numbers are as reported in the papers; see them for definitions, samples and caveats.

## Repository layout

| Path | Content |
|---|---|
| `paper_en.tex` | English paper (main text) |
| `paper_vi.tex` | Vietnamese version, with extended appendices |
| `preamble_common.tex` | Shared LaTeX preamble |
| `references.bib`, `references_pb.bib` | Bibliographies (Biber) |
| `figures/` | Figures used by the papers |
| `code/` | R and Python robustness checks and extra figures |
| `robust_out/` | Robustness tables (CSV) and plots |
| `docs/` | GitHub Pages site and compiled PDFs (`docs/pdf/`) |

## Build

XeLaTeX -> Biber -> XeLaTeX -> XeLaTeX, from the repository root:

```
xelatex paper_en && biber paper_en && xelatex paper_en && xelatex paper_en
```

Fonts: Times New Roman if installed, otherwise TeX Gyre Termes.

## Data and reproducibility

VHLSS microdata are licensed by the General Statistics Office and cannot be redistributed; no microdata, credentials or API tokens are in this repository. Scripts in `code/` need the licensed data and public climate inputs to rerun. The robustness script reads the data path from local configuration.

## Status

Working paper, not peer reviewed. Results may change.

## Citation

Truong, N. D. (2026). *Why Subsidised Aquaculture Insurance Does Not Reach Farmers: Climate Hazards, Trigger Design and Risk-Based Pricing in Vietnam's Mekong Delta.* Working paper. See `CITATION.cff`.

## Licence

Text, figures and PDFs: CC BY 4.0. Code: MIT. See `LICENSE`.