\* Problem: max_support *\

Maximize
 obj: - R + S

Subject To
 r_1: + alpha3 + alpha2 >= 1
 r_2: + alpha1 >= 1
 r_3: + 5 r1 + alpha3 <= 5
 r_4: + 5 r2 <= 5
 r_5: + 5 r3 + alpha4 + alpha3 <= 5
 r_6: + 5 r4 + alpha3 + alpha2 <= 5
 r_7: + 5 r5 + alpha1 <= 5
 r_8: - R + alpha4 + alpha2 + alpha1 = 0
 r_9: - R + alpha4 + alpha3 + alpha2 + alpha1 = 0
 r_10: + R <= 4
 r_11: - S + r5 + r4 + r3 + r2 + r1 = 0

Bounds
 0 <= alpha1 <= 1
 0 <= alpha2 <= 1
 0 <= alpha3 <= 1
 0 <= alpha4 <= 1
 0 <= r1 <= 1
 0 <= r2 <= 1
 0 <= r3 <= 1
 0 <= r4 <= 1
 0 <= r5 <= 1
 0 <= R <= 4
 0 <= S <= 5

Generals
 alpha1
 alpha2
 alpha3
 alpha4
 r1
 r2
 r3
 r4
 r5
 R
 S

End
