void ApplySettings()
{
	lhc_fill = 5277;

	cut1_a = -0.996; cut1_c = -0.797; cut1_si = 0.20;
	cut2_a = -0.972; cut2_c = +0.472; cut2_si = 0.20;
	cut3_a = -1.153; cut3_c = +0.087; cut3_si = 0.15;
	cut4_a = -1.147; cut4_c = +0.030; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(10.0, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.0, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 11.0);
}
