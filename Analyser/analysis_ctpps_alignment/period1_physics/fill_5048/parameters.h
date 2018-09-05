void ApplySettings()
{
	lhc_fill = 5048;

	cut1_a = -1.018; cut1_c = -0.775; cut1_si = 0.20;
	cut2_a = -0.988; cut2_c = +0.687; cut2_si = 0.20;
	cut3_a = -1.158; cut3_c = +0.315; cut3_si = 0.15;
	cut4_a = -1.143; cut4_c = +0.035; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(13.3, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(12.0, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(10.8, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(10.0, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(9.0, 12.5);
	yAlignmentRangeX["L_1_N"] = SelectionRange(9.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(6.0, 12.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
