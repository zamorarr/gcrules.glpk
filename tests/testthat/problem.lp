\* Problem: min_set_cover *\

Minimize
 obj: + R

Subject To
 r_1: + alpha3 + alpha2 >= 1
 r_2: + alpha1 >= 1
 r_3: - R + alpha4 + alpha2 + alpha1 = 0
 r_4: - R + alpha4 + alpha3 + alpha2 + alpha1 = 0

Bounds
 0 <= alpha1 <= 1
 0 <= alpha2 <= 1
 0 <= alpha3 <= 1
 0 <= alpha4 <= 1
 1 <= R <= 4

Generals
 alpha1
 alpha2
 alpha3
 alpha4
 R

End
