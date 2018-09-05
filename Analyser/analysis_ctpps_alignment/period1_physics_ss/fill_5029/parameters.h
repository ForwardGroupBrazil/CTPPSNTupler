void ApplySettings()
{
	lhc_fill = 5029;

	cut1_a = -1.014; cut1_c = -1.009; cut1_si = 0.20;
	cut2_a = -0.988; cut2_c = +0.664; cut2_si = 0.20;
	cut3_a = -1.163; cut3_c = +0.586; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.050; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.2, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.1, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.2, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(9.6, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(8.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.5, 12.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
