void ApplySettings()
{
	lhc_fill = 5045;

	cut1_a = -1.015; cut1_c = -0.842; cut1_si = 0.20;
	cut2_a = -0.985; cut2_c = +0.648; cut2_si = 0.20;
	cut3_a = -1.160; cut3_c = +0.345; cut3_si = 0.15;
	cut4_a = -1.143; cut4_c = +0.042; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.4, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.6, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.0, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(9.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.5, 12.5);
	yAlignmentRangeX["L_1_N"] = SelectionRange(8.5, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(6.0, 12.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(6.0, 12.0);
}
