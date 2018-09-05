void ApplySettings()
{
	lhc_fill = 5013;

	cut1_a = -0.999; cut1_c = -1.556; cut1_si = 0.20;
	cut2_a = -0.970; cut2_c = +0.434; cut2_si = 0.20;
	cut3_a = -1.163; cut3_c = +0.889; cut3_si = 0.15;
	cut4_a = -1.146; cut4_c = +0.038; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.5, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.0, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.0, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(7.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(4.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
