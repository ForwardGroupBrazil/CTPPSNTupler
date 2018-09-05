void ApplySettings()
{
	lhc_fill = 5267;

	cut1_a = -0.995; cut1_c = -0.807; cut1_si = 0.20;
	cut2_a = -0.969; cut2_c = +0.429; cut2_si = 0.20;
	cut3_a = -1.149; cut3_c = +0.076; cut3_si = 0.15;
	cut4_a = -1.136; cut4_c = +0.037; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(9.2, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(8.3, 15.7);
	selectionRangesX["R_1_N"] = SelectionRange(12.5, 16.0);
	selectionRangesX["R_1_F"] = SelectionRange(12.0, 15.1);

	yAlignmentRangeX["L_1_F"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["L_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(4.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
