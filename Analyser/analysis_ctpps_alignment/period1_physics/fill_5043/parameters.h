void ApplySettings()
{
	lhc_fill = 5043;

	cut1_a = -1.013; cut1_c = -0.898; cut1_si = 0.20;
	cut2_a = -0.982; cut2_c = +0.580; cut2_si = 0.20;
	cut3_a = -1.165; cut3_c = +0.398; cut3_si = 0.15;
	cut4_a = -1.142; cut4_c = +0.033; cut4_si = 0.15;

	selectionRangesX["L_1_F"] = SelectionRange(12.0, 16.5);
	selectionRangesX["L_1_N"] = SelectionRange(11.2, 14.5);
	selectionRangesX["R_1_N"] = SelectionRange(9.5, 15.5);
	selectionRangesX["R_1_F"] = SelectionRange(8.8, 14.5);

	yAlignmentRangeX["L_1_F"] = SelectionRange(8.0, 12.5);
	yAlignmentRangeX["L_1_N"] = SelectionRange(8.0, 12.0);
	yAlignmentRangeX["R_1_N"] = SelectionRange(5.0, 11.0);
	yAlignmentRangeX["R_1_F"] = SelectionRange(4.0, 11.0);
}
