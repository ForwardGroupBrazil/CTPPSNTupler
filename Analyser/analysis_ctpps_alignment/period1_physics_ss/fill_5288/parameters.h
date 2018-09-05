void ApplySettings()
{
	lhc_fill = 5288;

	cut1_a = -1.008; cut1_c = -0.598; cut1_si = 0.20;
	cut2_a = -0.992; cut2_c = +0.730; cut2_si = 0.20;
	cut3_a = -1.143; cut3_c = +0.157; cut3_si = 0.15;
	cut4_a = -1.130; cut4_c = +0.063; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(10.1, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.7, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(7.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 11.0);
}
