int cholmod_suite_setup(void);
int cholmod_suite_teardown(void);
void cholmod_test_setup(void);
void cholmod_test_teardown(void);

void test_mat_vec(void);
void test_mat_tpose_vec(void);
void test_mat_inf_norm_cols(void);
void test_mat_inf_norm_rows(void);
void test_ldlchol(void);
