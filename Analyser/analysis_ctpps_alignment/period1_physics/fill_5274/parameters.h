void ApplySettings()
{
	lhc_fill = 5274;

	cut1_a = -0.994; cut1_c = -0.799; cut1_si = 0.20;
	cut2_a = -0.970; cut2_c = +0.436; cut2_si = 0.20;
	cut3_a = -1.158; cut3_c = +0.088; cut3_si = 0.15;
	cut4_a = -1.148; cut4_c = +0.024; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.5, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.5, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(4.5, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
