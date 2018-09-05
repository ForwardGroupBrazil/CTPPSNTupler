void ApplySettings()
{
	lhc_fill = 5276;

	cut1_a = -0.995; cut1_c = -0.800; cut1_si = 0.20;
	cut2_a = -0.971; cut2_c = +0.440; cut2_si = 0.20;
	cut3_a = -1.150; cut3_c = +0.094; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.030; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.5, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.8, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(4.5, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
