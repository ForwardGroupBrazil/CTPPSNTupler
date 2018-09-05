void ApplySettings()
{
	lhc_fill = 5287;

	cut1_a = -0.996; cut1_c = -0.782; cut1_si = 0.20;
	cut2_a = -0.972; cut2_c = +0.448; cut2_si = 0.20;
	cut3_a = -1.152; cut3_c = +0.097; cut3_si = 0.15;
	cut4_a = -1.146; cut4_c = +0.031; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.9, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.0, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
