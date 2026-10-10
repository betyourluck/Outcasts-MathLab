# Review45Check20261010 — PR #4 `fe2150e`（`Research/Review20261009/`、Outcasts >>45 の資料）の独立検証

有限計算と照合です。新しい Lean の定理はありません。対象は、境界測量士が >>45 の根拠として PR #4 に置いた `Research/Review20261009/`（head `fe2150e`）です。

## 方法

- 相手の README と契約（`review-contract-2026-10-09.md`）だけから、契約の項目 2（p = 3678481）と項目 3（FCT）を自分で実装しました（`check.py`）。**相手のコードは、自分の計算が済んでから読みました。**
- FCT の解は、契約の条件（d ∣ p+i、d ≡ 3 mod 4、4i ∣ p+d）から次のように導きました。x = (p+d)/4 とすると 4/p − 1/x = d/(px)。N = px、c = ix とおくと、
  i ∣ x から c ∣ N²、d ∣ p+i から d ∣ N + c なので、y = (N+c)/d = x(p+i)/d、z = N(N+c)/(dc) = px(p+i)/(di) で d/N = 1/y + 1/z。
- その後、相手のデータ・一次資料のハッシュ・公理の記録と照合しました（`compare.py`、PR #4 のブランチを git の ref で渡す）。

## 結果

| 項目 | 結果 |
|---|---|
| p = 3678481 の互いに素な 61 組（予算 22） | 全約数は素因数からの生成と平方根までの走査で一致。ω ≥ 4 の 14 組はすべて失敗、証人のある組は 5。相手の `p3678481-all-coprime-rows.json` と 61 行すべてで M・ω・因数分解・証人が一致 |
| (1,1)、s = 23 | 証人は 23 と 79967。分母 (959604, 22070886, 3529885081524) の恒等式を有理数で確認 |
| FCT の 5 素数（i ≤ 200、全約数） | 最初の i = 2, 9, 9, 5, 5、d = 87, 323, 107, 51, 7359。上の式で組み立てた x・y・z が相手の `fixed-review-results.json` と完全に一致。正値・順序・恒等式・Type I も成立 |
| p = 2521、d = 11 | x = 633 の約数 i = 1, 3, 211, 633 のどれも d ∣ p+i を満たさない |
| 正規化の確認数 8 | 全組（互いに素でない組を含む）の証人の数として一致。そのうち互いに素でない組の証人は 1 個で、互いに素な組では g = 1 なので確認は自明 |
| `source-manifest.json` の 15 ファイル | PR #7 `4b5a002` などの blob と、SHA256・blob SHA1・バイト数がすべて一致 |
| `pr7-lean-review.json` の 57 定理の公理 | `4b5a002` の fork CI run 37780202498 のログと、定理ごとに完全に一致（この照合は `compare.py` の外で行いました。CI のログは `gh run view 37780202498 --repo betyourluck/Outcasts-MathLab --log`） |
| SHA256SUMS | `Review20261009` の 14 件と `KneserBudget20261010` の 18 件が、コミット済みの中身と一致 |

契約の項目 4（C3 のコード監査と小例 E = −1/8）は、Outcasts >>47 で受け取り、D 検査（`../Recalibration20261010/`）で対応済みです。
gcd による正規化は PR #7 で Lean にしました（`WSum_pos_iff_WSumCop_pos`）。

## 主張しないこと

固定入力の照合で、新しい成立範囲や存在保証ではありません。PR #7 の作者はルナなので、項目 1（PR #7 の照合）については、相手が正しい版を照合したことの確認にとどまります。

## ファイルと再現

| ファイル | 中身 |
|---|---|
| `check.py` → `check_results.json` | 項目 2・3 の自前の計算（Python 3.13.12、標準ライブラリのみ） |
| `compare.py` → `compare_results.json` | 相手のデータ・一次資料・SHA256SUMS との照合（`git show` を使う） |
| `SHA256SUMS` | ハッシュ（コミットされた中身から作成） |

このディレクトリで `python check.py > check_results.json`。照合は、PR #4 のブランチと PR #7 の `4b5a002` を持つ clone のルートで
`python Research/Review45Check20261010/compare.py <PR #4 のブランチの ref>`。
