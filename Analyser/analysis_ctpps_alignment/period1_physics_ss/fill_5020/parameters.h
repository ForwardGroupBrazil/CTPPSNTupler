void ApplySettings()
{
	lhc_fill = 5020;

	cut1_a = -1.009; cut1_c = -1.336; cut1_si = 0.20;
	cut2_a = -0.979; cut2_c = +0.566; cut2_si = 0.20;
	cut3_a = -1.165; cut3_c = +0.778; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.032; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(11.0, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(9.3, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(8.9, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(8.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(6.5, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(6.5, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 10.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(5.0, 10.0);
}
