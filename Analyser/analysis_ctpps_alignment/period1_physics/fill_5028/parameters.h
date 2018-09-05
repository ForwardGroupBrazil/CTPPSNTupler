void ApplySettings()
{
	lhc_fill = 5028;

	cut1_a = -1.014; cut1_c = -1.053; cut1_si = 0.20;
	cut2_a = -0.985; cut2_c = +0.640; cut2_si = 0.20;
	cut3_a = -1.161; cut3_c = +0.589; cut3_si = 0.15;
	cut4_a = -1.143; cut4_c = +0.037; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.1, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(10.7, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.1, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(9.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(7.5, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
