void ApplySettings()
{
	lhc_fill = 5045;

	cut1_a = -1.021; cut1_c = -0.766; cut1_si = 0.20;
	cut2_a = -0.988; cut2_c = +0.688; cut2_si = 0.20;
	cut3_a = -1.159; cut3_c = +0.347; cut3_si = 0.15;
	cut4_a = -1.144; cut4_c = +0.038; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.4, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.6, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.0, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(9.5, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(9.0, 12.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(9.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(7.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(6.0, 11.0);
}
