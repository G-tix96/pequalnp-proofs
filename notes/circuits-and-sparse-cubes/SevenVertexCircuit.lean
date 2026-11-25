import Mathlib
import PEqualNP.Defs

open PEqualNP.CircuitLowerBounds

private def witness : Circuit := ⟨[
  ⟨.or, [0, 1]⟩, ⟨.or, [2, 3]⟩, ⟨.or, [5, 6]⟩, ⟨.and, [7, 4]⟩,
  ⟨.or, [7, 8]⟩, ⟨.or, [11, 4]⟩, ⟨.and, [9, 12]⟩, ⟨.or, [10, 13]⟩], 14⟩

example : witness.gates.length = 8 := rfl

example : ∀ x : Fin 7 → Bool, witness.eval (List.ofFn x) = decide (∃ i < ([0, 4, 1, 4, 0, 5, 1, 5, 2, 5, 3, 5, 4, 5, 0, 6, 1, 6, 2, 6, 3, 6, 4, 6] : List ℕ).length / 2, (List.ofFn x).getD (([0, 4, 1, 4, 0, 5, 1, 5, 2, 5, 3, 5, 4, 5, 0, 6, 1, 6, 2, 6, 3, 6, 4, 6] : List ℕ).getD (2 * i) 0) false = true ∧ (List.ofFn x).getD (([0, 4, 1, 4, 0, 5, 1, 5, 2, 5, 3, 5, 4, 5, 0, 6, 1, 6, 2, 6, 3, 6, 4, 6] : List ℕ).getD (2 * i + 1) 0) false = true) := by
  decide +kernel
