# -*- coding: utf-8 -*-
"""paper_extra_figures.py (18/10/2026): hai hình bổ sung cho paper_vi.tex.
jv13_sensitivity.pdf : T28_sensitivity_supply.csv -> mức tăng tỷ lệ tham gia của hộ không nghèo khi chuyển từ phí chung sang phí theo ô.
jv14_trigger_lift.pdf: Bảng tab:auc (giá trị đã công bố trong bài) -> hiệu quả ngưỡng phân vị 80 so với ngưỡng ngẫu nhiên.
Chạy: python -I paper_extra_figures.py <T28_sensitivity_supply.csv> <thư mục đầu ra>"""
import sys, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in ["/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf",
          "/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf"]:
    if os.path.exists(f): font_manager.fontManager.addfont(f)
plt.rcParams.update({"font.family": "TeX Gyre Termes", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "legend.frameon": False})
csv, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(csv)
# hàng 'flat' không ghi s_oth; suy ra từ thứ tự: tại mỗi (load) hai hàng flat, hàng có tham gia thấp hơn là s=0,2
rows = []
for load, g in d.groupby("load"):
    fl = g[g.scheme == "flat"].sort_values("takeup_nonpoor")
    flat = {0.2: fl.iloc[0].takeup_nonpoor, 0.5: fl.iloc[1].takeup_nonpoor}
    for _, r in g[g.scheme != "flat"].iterrows():
        rows.append(dict(load=load, rho=r.rho, s=r.s_oth, gain=100 * (r.takeup_nonpoor - flat[r.s_oth]),
                         flat=100 * flat[r.s_oth], map=100 * r.takeup_nonpoor))
t = pd.DataFrame(rows)
t.to_csv(os.path.join(out, "jv13_sensitivity_data.csv"), index=False)
C = {0.5: "#9DB4D8", 0.8: "#3B6BB0", 0.95: "#00196E"}
fig, ax = plt.subplots(1, 2, figsize=(6.3, 2.7), sharey=True)
for a, s, ttl in zip(ax, [0.2, 0.5], ["(a) Hỗ trợ hộ không nghèo 20%", "(b) Hỗ trợ hộ không nghèo 50%"]):
    for rho in [0.5, 0.8, 0.95]:
        x = t[(t.s == s) & (t.rho == rho)].sort_values("load")
        a.plot(100 * x.load, x.gain, marker="o", ms=4, lw=1.4, color=C[rho], label=f"$\\rho$ = {rho}".replace(".", ","))
        for _, r in x.iterrows():
            a.annotate(f"{r.gain:.1f}".replace(".", ","), (100 * r.load, r.gain), textcoords="offset points",
                       xytext=(0, 5), ha="center", fontsize=7, color=C[rho])
    a.set_title(ttl, fontsize=9, loc="left")
    a.set_xlabel("Phụ phí quản lý (%)")
    a.set_xticks([15, 25, 40]); a.set_ylim(0, 32); a.grid(axis="y", lw=0.3, color="#cccccc")
ax[0].set_ylabel("Tăng tỷ lệ tham gia của hộ không\nnghèo so với phí chung (điểm %)")
ax[0].legend(loc="upper right", fontsize=8, title="Độ chính xác bản đồ", title_fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(out, "jv13_sensitivity.pdf"), bbox_inches="tight"); plt.close(fig)

# ---- hình lift: số liệu từ bảng tab:auc của bài -------------------------------------------
L = [("Hạn: ngày Tmax ≥35°C", .19, .96), ("Hạn: SPEI-3 nhỏ nhất", .37, .92), ("Hạn: mưa mùa khô", .32, .94),
     ("Hạn: khô–nóng (PC1)", .21, .96), ("Lũ: ngày mưa ≥50 mm", .50, .80), ("Lũ: mưa ngày lớn nhất", .57, .77),
     ("Lũ: dị thường mưa 12 tháng", .35, .87), ("Lũ: chỉ số ướt (PC1)", .46, .82), ("Bão: chỉ số phá hủy", .21, .90)]
fig, ax = plt.subplots(1, 2, figsize=(6.3, 3.0), sharey=True, gridspec_kw=dict(width_ratios=[1, 1]))
y = list(range(len(L)))[::-1]
col = ["#7F7F7F"] * 4 + ["#00196E"] * 4 + ["#9DB4D8"]
ax[0].barh(y, [l[1] for l in L], color=col, height=0.62)
ax[0].axvline(0.20, color="#C0392B", lw=1, ls="--")
for yy, l in zip(y, L): ax[0].text(l[1] + 0.012, yy, f"{l[1]/0.20:.2f}".replace(".", ",") + "×", va="center", fontsize=7.5)
ax[0].set_yticks(y); ax[0].set_yticklabels([l[0] for l in L], fontsize=8)
ax[0].set_xlim(0, 0.75); ax[0].set_xlabel("(a) Tỷ lệ xã báo sự kiện được chi trả"); 
ax[0].text(0.205, y[0] + 0.62, "ngẫu nhiên = 0,20", color="#C0392B", fontsize=7, va="bottom")
ax[1].barh(y, [l[2] for l in L], color=col, height=0.62)
for yy, l in zip(y, L): ax[1].text(l[2] + 0.012, yy, f"{100*l[2]:.0f}%", va="center", fontsize=7.5)
ax[1].set_xlim(0, 1.15); ax[1].set_xlabel("(b) Tỷ lệ lượt chi trả rơi vào xã không báo sự kiện")
ax[1].tick_params(left=False)
fig.tight_layout(); fig.savefig(os.path.join(out, "jv14_trigger_lift.pdf"), bbox_inches="tight"); plt.close(fig)
print(t.gain.min(), t.gain.max())
