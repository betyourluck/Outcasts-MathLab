"""Supplement (not in the contract, written after D1-D4 were computed): first-order size of the
same-sample refit effect on E[T2] under the null of D1.

Within a calibration cell c, P''_i = P_i * O*/E_c with E_c = sum_c P_k.  Writing eps_k = f*_k - P'_k
(independent, variance v_k = P'_k(1 - P'_k)) and ignoring clipping,
    d_i = f*_i - P''_i = eps_i - (P_i / E_c) * sum_{k in c} eps_k,
so for i != j in the same cell
    E[d_i d_j] = -(P_j v_i + P_i v_j) / E_c + P_i P_j (sum_{k in c} v_k) / E_c^2,
and 0 for different cells.  E[T2] is the sum over pairs at the same p.  (For two rows with
P = 1/2 in one cell this gives -1/8, the example in Outcasts >>45.)"""
import gzip, json, sys
from collections import defaultdict
from fractions import Fraction

import recalibration_test as rt

def bias(rows, cell, Pkey):
    E = defaultdict(float); V = defaultdict(float)
    for r in rows:
        if r["method"] == "deterministic": continue
        c = cell(r); E[c] += r["P"]; V[c] += r[Pkey] * (1 - r[Pkey])
    groups = defaultdict(list)
    for r in rows:
        if r["method"] != "deterministic": groups[(r["p"], cell(r))].append(r)
    total = 0.0
    for (p, c), g in groups.items():
        if E[c] == 0: continue  # the "P = 0" cell: P'' = 0 and f* = 0, no contribution
        for x in range(len(g)):
            for y in range(x + 1, len(g)):
                i, j = g[x], g[y]
                vi, vj = i[Pkey] * (1 - i[Pkey]), j[Pkey] * (1 - j[Pkey])
                total += -(j["P"] * vi + i["P"] * vj) / E[c] + i["P"] * j["P"] * V[c] / E[c] ** 2
    return total

def toy():
    rows = [{"p": 1, "P": Fraction(1, 2), "Pp": Fraction(1, 2), "method": "exact"} for _ in range(2)]
    E = sum(r["P"] for r in rows); V = sum(r["Pp"] * (1 - r["Pp"]) for r in rows)
    i, j = rows
    return str(-(j["P"] * i["Pp"] * (1 - i["Pp"]) + i["P"] * j["Pp"] * (1 - j["Pp"])) / E + i["P"] * j["P"] * V / E ** 2)

def main():
    conf = json.loads(gzip.open(rt.jt_path("confirmatory_rows.json.gz")).read())
    rows = rt.st.model_rows(conf, 23, 10000, 20261008)
    rt.apply(rows, rt.ratios(rows, rt.cell_C3), rt.cell_C3, "P_recal")
    out = {"toy_formula": toy(),
           "E_T2_shift_method_x_Pbin": round(bias(rows, rt.cell_C3, "P_recal"), 3)}
    rt.apply(rows, rt.ratios(rows, rt.cell_box), rt.cell_box, "P_box")
    out["E_T2_shift_method_x_Pbin_x_box"] = round(bias(rows, rt.cell_box, "P_box"), 3)
    json.dump(out, sys.stdout, indent=1)

if __name__ == "__main__":
    main()
