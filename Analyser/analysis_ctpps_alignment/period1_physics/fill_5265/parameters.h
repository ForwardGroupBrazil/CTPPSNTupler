void ApplySettings()
{
	lhc_fill = 5265;

	cut1_a = -0.992; cut1_c = -0.825; cut1_si = 0.20;
	cut2_a = -0.971; cut2_c = +0.451; cut2_si = 0.20;
	cut3_a = -1.142; cut3_c = +0.100; cut3_si = 0.15;
	cut4_a = -1.140; cut4_c = +0.035; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(8.7, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.3, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(4.0, 10.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
