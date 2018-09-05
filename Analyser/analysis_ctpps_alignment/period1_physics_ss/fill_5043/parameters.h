void ApplySettings()
{
	lhc_fill = 5043;

	cut1_a = -1.007; cut1_c = -0.972; cut1_si = 0.20;
	cut2_a = -0.980; cut2_c = +0.557; cut2_si = 0.20;
	cut3_a = -1.165; cut3_c = +0.410; cut3_si = 0.15;
	cut4_a = -1.145; cut4_c = +0.032; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.0, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.2, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(9.5, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(8.8, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(7.5, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(7.5, 11.5);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 11.0);
}
