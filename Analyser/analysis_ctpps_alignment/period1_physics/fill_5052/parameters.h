void ApplySettings()
{
	lhc_fill = 5052;

	cut1_a = -1.019; cut1_c = -0.710; cut1_si = 0.20;
	cut2_a = -0.989; cut2_c = +0.678; cut2_si = 0.20;
	cut3_a = -1.156; cut3_c = +0.298; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.039; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(13.3, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(12.0, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.8, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(10.1, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(9.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(9.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
