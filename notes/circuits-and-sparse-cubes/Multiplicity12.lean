import Mathlib

/- A local lower bound for the family in L-gct-237. -/
namespace SparseCubeExperiment

open Polynomial

noncomputable section

private def tri (a b c : ℝ) (i j k : ℕ) : Polynomial ℝ :=
  C a * X ^ i + C b * X ^ j + C c * X ^ k

private def A (r : ℝ) := tri (143 * r) (-286 * r) (143 * r) 1 8 15
private def B := tri 98 (-637) 539 0 11 13
private def D := tri 539 (-637) 98 3 5 16

private lemma sparse (a b c : ℝ) (i j k : ℕ) : (tri a b c i j k).support.card ≤ 3 := by
  apply (Finset.card_le_card (support_trinomial_subset i j k a b c)).trans
  exact Finset.card_le_three

-- The quotient by (X - 1)^12 is reciprocal of degree 36.
private def Q : Polynomial ℝ :=
  941192 * (1 + X^36) + 11294304 * (X + X^35) + 73412976 * (X^2 + X^34) +
  336745474 * (X^3 + X^33) + 1214546112 * (X^4 + X^32) + 3654950364 * (X^5 + X^31) +
  9519369496 * (X^6 + X^30) + 21969409098 * (X^7 + X^29) + 45591301392 * (X^8 + X^28) +
  85859227475 * (X^9 + X^27) + 147767730648 * (X^10 + X^26) + 233973758811 * (X^11 + X^25) +
  343119030904 * (X^12 + X^24) + 468976601226 * (X^13 + X^23) + 600789791196 * (X^14 + X^22) +
  724781584254 * (X^15 + X^21) + 826469068152 * (X^16 + X^20) +
  893220184155 * (X^17 + X^19) + 916477221244 * X^18

private lemma factorization : B^3 + D^3 - 2 * (A 1)^3 = Q * (X - C 1)^12 := by
  norm_num [A, B, D, tri, Q, C_ofNat]
  ring

theorem exists_order_twelve :
    ∃ (s : Fin 3 → ℝ) (g : Fin 3 → Polynomial ℝ),
      (∀ i, s i = 1 ∨ s i = -1) ∧ (∀ i, (g i).support.card ≤ 3) ∧
      (∑ i, C (s i) * g i ^ 3) ≠ 0 ∧
      rootMultiplicity 1 (∑ i, C (s i) * g i ^ 3) = 12 := by
  let r := (2 : ℝ) ^ ((3 : ℝ)⁻¹)
  have hr : r ^ 3 = 2 := Real.rpow_inv_natCast_pow (by norm_num) (by decide : 3 ≠ 0)
  have ha : A r = C r * A 1 := by
    simp [A, tri]
    ring
  let s : Fin 3 → ℝ := ![-1, 1, 1]
  let g : Fin 3 → Polynomial ℝ := ![A r, B, D]
  refine ⟨s, g, ?_, ?_, ?_⟩
  · intro i
    fin_cases i <;> simp [s]
  · intro i
    fin_cases i <;> apply sparse
  · have hsum : (∑ i, C (s i) * g i ^ 3) = Q * (X - C 1)^12 := by
      rw [← factorization]
      simp [s, g, Fin.sum_univ_succ, ha, mul_pow, ← map_pow, hr, C_ofNat]
      ring
    have hval : Q.eval 1 ≠ 0 := by norm_num [Q, eval_add, eval_mul, eval_pow]
    have hQ : Q ≠ 0 := by intro h; exact hval (by simp [h])
    rw [hsum, rootMultiplicity_mul_X_sub_C_pow hQ, rootMultiplicity_eq_zero hval]
    exact ⟨mul_ne_zero hQ (pow_ne_zero _ (X_sub_C_ne_zero 1)), rfl⟩

end
end SparseCubeExperiment
