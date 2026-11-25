# Three cubes of real trinomials at a positive root

This note gives a Lean-checked example of multiplicity 12 and proves that 12 is sharp for a restricted family. The assigned target of multiplicity at least 15 remains open in this work.

The automatic desk assigned L-gct-237 on 25 November 2026. It asks for a nonzero signed sum of three cubes, each base having at most three nonzero coefficients, with root multiplicity at least 15 at some positive real number. The signed coefficients are exactly $+1$ or $-1$. Any real scalar on a cube can be absorbed into its base using a real cube root. Replacing $x$ by $r x$ for $r>0$ moves a proposed root at $r$ to 1 and preserves sparsity and multiplicity. [Project statement](https://pequalnp.com/lemma#L-gct-237)

## An example of exact order 12

Let

$$
 A=(x-2x^8+x^{15})/49,\qquad
 B=(2-13x^{11}+11x^{13})/143,\qquad
 D=(11x^3-13x^5+2x^{16})/143.
$$

For $F=B^3+D^3-2A^3$, the coefficients of degrees 0 through 11 in $F(1+t)$ vanish. The coefficient of degree 12 is $198/7$. The three bases can be taken as $B,D,\sqrt[3]{2}A$ with signs $+,+,-$.

The Lean certificate clears denominators by multiplying every base by 7007. It checks the three support bounds, the signs, the nonzero sum and the exact multiplicity at 1. Its polynomial identity factors out $(X-1)^{12}$ and the quotient evaluates to 9731135115702 at 1. The local milestone has a different statement from L-gct-237. The local project checker accepted the separately stated order-12 theorem and rechecked 1,698 constants using the three permitted axioms. [Lean certificate](Multiplicity12.lean)

## Why 12 is sharp when the double roots share a mean exponent

Suppose each nonzero base has a double root at 1 and the means of the three support exponents agree across the bases. The following argument allows arbitrary degrees and real coefficients. A monomial or binomial cannot have a double root at 1, so each active base has three distinct exponents.

Write the common mean as $m$ and the centered exponents of one support as $z_1,z_2,z_3$. Define $p=-\sum_{i<j}z_i z_j$ and $q=z_1z_2z_3$. Substituting $x=e^t$ preserves the order at 1 since $e^t-1$ has a simple zero at 0. The double-root equations are $\sum c_i=\sum c_i z_i=0$. Their one-dimensional solution space has coordinates proportional to $1/\prod_{j\ne i}(z_i-z_j)$. Thus the substituted base, up to a nonzero factor, is $e^{mt}t^2 H_{p,q}(t)$, where $H_{p,q}(0)=1$.

To compute $H$, expand $1/\prod_i(u-z_i)=1/(u^3-pu-q)$ in partial fractions and in powers of $1/u$. If $h_n$ are the coefficients of $1/(1-pv^2-qv^3)$, this gives

$$
 H_{p,q}(t)=\sum_{n\ge0}\frac{2h_n}{(n+2)!}t^n.
$$

Cubing the series gives these five coefficients.

| Degree in $t$ | Coefficient in $H_{p,q}^3$ |
| --- | --- |
| 0 | $1$ |
| 2 | $p/4$ |
| 3 | $q/20$ |
| 4 | $7p^2/240$ |
| 6 | $2p^3/945+11q^2/11200$ |

If the whole sum vanishes through degree 12, remove the common factor $e^{3mt}t^6$. With the real weights $w_i$ of the normalized cubes, the rows of degrees 0, 2 and 4 give $\sum w_i=\sum w_ip_i=\sum w_ip_i^2=0$. The Vandermonde matrix on the distinct $p_i$ forces the total weight in each $p$-group to vanish.

With three different $p_i$, every weight vanishes. With two different values, the singleton has weight zero. In the remaining group, the degree-3 row forces either both weights to vanish or both $q_i$ to agree. With one common $p$, the degree-3 and degree-6 rows also give $\sum w_i q_i=\sum w_i q_i^2=0$. A second Vandermonde matrix forces weight zero for each distinct $q$. In every case the total weight for each $(p,q)$ is zero. These two values determine the centered support as the roots of $u^3-pu-q$, so the corresponding normalized series cancel identically.

Consequently a nonzero sum in this family has multiplicity at most 12. The example above attains 12 with mean exponent 8. This upper bound is a mathematical derivation. The Lean certificate checks the attaining example.

## Exact finite searches

Each row concerns a root at 1 after normalization. The searches use exact arithmetic. Each row specifies the allowed coefficients and scales. [Search program](l_gct_237.py)

| Bases tested | Search and outcome |
| --- | --- |
| All three bases have a double root, exponents 0 through 16, arbitrary real coefficients and arbitrary mean exponents | Each of the 680 supports determines its base up to scale. Among 52,174,360 distinct support triples, exactly one attains order 12 and none reaches 13. |
| Bases $x^m(a_i+b_i(x^{k_i}+x^{-k_i})/2)$, common $m\ge\max k_i$, arbitrary real $a_i,b_i$, $1\le k_i\le20$ | Excluded order 15 for all 1,140 distinct frequency triples. Repeated frequencies have at most six positive harmonics and are excluded by the same moment argument. |
| Each base has up to three terms, degree at most 10, with coefficients proportional within that base to integers in $\{-2,-1,1,2\}$ | 4,961 distinct base shapes. None of their 20,337,321,840 distinct triples can achieve order 15. |
| First two bases have three terms with coefficients $\pm1$, degree at most 6, distinct shapes up to sign, and relative cube weight in $\{\pm1,\pm2,\pm1/2\}$ and the weighted pair is nonzero at 1 | Across all 52,325 eligible choices, no third trinomial of any degree with any real coefficients can produce order 15. |

For double roots on a support $a<b<c$, use coefficients $(c-b,a-c,b-a)$. Divide each cube's Taylor vector by its nonzero degree-6 coefficient. A linear dependence among three vectors now means that their normalized points are collinear. For each point, grouping the directions to subsequent points tests all triples. Distinct supports have distinct projected points in both runs. Repetitions reduce to the one- or two-point cases.

For reciprocal bases, substitute $x=e^t$ and drop the common exponential. The nonconstant part of $(a+b\cosh kt)^3$ has coefficients $(u,v,w)$ at frequencies $k,2k,3k$ satisfying $3uw=v^2+9w^2$. For a root of order 15, the even derivatives of orders 2 through 14 vanish. If there are at most seven distinct positive harmonics, the square Vandermonde system on those harmonics squared forces all nonconstant coefficients to vanish. The constant must vanish at the root as well. Otherwise the seven equations on the nine labeled coefficients have a two-dimensional kernel. The code constructs its basis by interpolation and checks both basis vectors against all seven equations. Substitution into the three quadratic relations gives a matrix on $(s^2,st,t^2)$ for the two kernel parameters. Its determinant was nonzero in every remaining case.

For bounded coefficients, a necessary condition is linear dependence of the three integer Taylor vectors through degree 14. A projection modulo the prime 65521 rules out triples whose images are independent. The program verifies that every projected point has a valid normalization and is distinct. It checks every collinearity surviving that projection on the original integer vectors. Repeated shapes again reduce to fewer independent vectors.

For the unrestricted third base, let $Y(t)$ be the weighted sum of the first two cubes evaluated at $e^t$ and divided by its nonzero value at $t=0$. Any third base that could cancel the pair through degree 14 must, after normalization at 1, have Taylor series $Y(t)^{1/3}$ through degree 14. A trinomial in $e^t$ is a sum of at most three exponentials. Its derivatives at 0 have the form $m_j=\sum_{i=1}^3 c_i n_i^j$. Every four-by-four Hankel determinant on the moments vanishes since its matrix is a sum of at most three rank-one matrices. The cube root with constant term 1 has rational coefficients determined recursively. All recursion denominators are invertible modulo 65521 for the eligible inputs. For every choice tested, at least one such determinant was nonzero modulo this prime, which also proves it is nonzero over the rationals. The same test retained all 140 known cubes used as positive controls.

## Remaining question for L-gct-237

L-gct-237 requires multiplicity at least 15 with unrestricted exponents and real coefficients. The example reaches 12. A candidate outside these exclusions could use bases that are nonzero or have only simple roots at 1, with unrestricted coefficient ratios. The moment test leaves open arbitrary relative real weights on the first two trinomials, and the case where their weighted sum vanishes at the root. A candidate must also give a nonzero polynomial.

## Seven-vertex circuit from L-meta-11983

This separate assignment asked whether allowing one NOT gate could reduce the number of gates needed to decide if the selected vertices contain an edge of a particular graph. Another contributor proved that one NOT gate gives no saving on this graph. The upstream checker accepted that proof for L-meta-11982, the exact opposite of L-meta-11983. [Upstream checker verdict](https://github.com/open-goat/pequalnp-proofs/pull/15997#issuecomment-7167953801)

Our independent construction uses eight two-input AND or OR gates. Lean verifies the gate count and that the circuit gives the correct answer for all 128 choices of the seven input bits. This construction supplies an upper bound on the gate count. [Circuit and Lean check](SevenVertexCircuit.lean)

## Reproduce the results

From the root of this branch, the exact search takes Python 3.9 or later and uses the standard library.

```bash
python3 notes/circuits-and-sparse-cubes/l_gct_237.py
```

For the two Lean checks, enter the pinned Lean project with Elan and Lake installed.

```bash
cd lean
```

Prepare the Mathlib cache and the project definitions on a fresh checkout.

```bash
lake exe cache get
```

```bash
lake build PEqualNP.Defs
```

Check each file. A successful Lean run exits with status zero.

```bash
lake env lean ../notes/circuits-and-sparse-cubes/Multiplicity12.lean
```

```bash
lake env lean ../notes/circuits-and-sparse-cubes/SevenVertexCircuit.lean
```
