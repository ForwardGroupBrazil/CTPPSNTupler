void ApplySettings()
{
	lhc_fill = 5030;

	cut1_a = -1.020; cut1_c = -0.915; cut1_si = 0.20;
	cut2_a = -0.990; cut2_c = +0.714; cut2_si = 0.20;
	cut3_a = -1.159; cut3_c = +0.500; cut3_si = 0.15;
	cut4_a = -1.143; cut4_c = +0.038; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.7, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.6, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.7, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(10.0, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.5, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(8.5, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
